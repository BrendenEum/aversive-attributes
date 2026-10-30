# Parameter recovery for the main Bayesian attribute-level model (variant
# cond_lam: one prior per range condition, shared loss aversion, per-subject
# noise partially pooled; see code/fitting/fit_bayesian_attributes.R).
#
# From the project root:
#
#   caffeinate -i nohup Rscript code/parameter_recovery/recover_bayesian_attributes.R > results/logs/recovery.log 2>&1 &
#
# For each of N_DATASETS datasets (20; each takes a few minutes): draw true
# parameters, simulate one choice per real trial (real stimuli, real dwell
# times), refit from scratch with the profile-likelihood procedure in
# code/fitting/fit_bayesian_attributes.R, and store true and recovered values.
# Trials with an unhovered box are dropped, as in the real fits. Results are
# checkpointed after every dataset to
# results/cache/bayesian_model/parameter_recovery.rds (re-running the script
# resumes from there, so with a complete cache it only redraws the figure), and
# the figure (plot_parameter_recovery() in code/utilities/) goes to
# paper/figures/parameter_recovery_bayesian.pdf.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source(here::here("code/utilities/getAllUtilities.R"))          # paths, .myPlot, build_pas(), build_boxlong()
source(here::here("code/simulation/simulate_bayesian_attributes.R"))
source(here::here("code/fitting/fit_bayesian_attributes.R"))

N_DATASETS <- 20
SEED       <- 2026
# Coarser profile grid than the real fit, to keep a dataset to ~1.5 minutes.
Q_GRID     <- c(-0.1, 0.1, 0.3, 0.5)

cache_dir  <- file.path(paths$root, paths$results$cache, "bayesian_model")
fig_dir    <- file.path(paths$root, paths$paper_figures)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
cache_file <- file.path(cache_dir, "parameter_recovery.rds")

# ---------------------------------------------------------------------------
# Real trial structure
# ---------------------------------------------------------------------------
pas_raw <- read_csv(file.path(paths$root, paths$data$processed, "main_data.csv"), show_col_types = FALSE)
pas     <- build_pas(pas_raw)
boxlong <- build_boxlong(pas)

variant  <- bayesian_variants()$cond_lam      # per-condition prior, shared loss aversion, partial pooling
prep     <- prepare_boxes(boxlong, variant)
subjects <- prep$subjects                      # Study, Subject, WN
conds    <- names(prep$cond)                   # "Narrow", "Wide"
# Sign of the true value difference per trial, to report simulated accuracy.
true_diff <- prep$all$trials %>% left_join(pas %>% select(Subject, Trial, Diff), by = c("Subject", "Trial"))

# ---------------------------------------------------------------------------
# True parameters, all drawn up front so that resuming reproduces the same draws
# ---------------------------------------------------------------------------
set.seed(SEED)
true_params <- list()
for (d in seq_len(N_DATASETS)) {
  true_params[[d]] <- list(
    sigma = setNames(exp(runif(nrow(subjects), log(2), log(25))), subjects$Subject),    # noise per subject, points per second (fitted values: median ~7, 5-95% ~2.5-17)
    s0    = setNames(rep(1, length(conds)), conds),                                     # prior SD fixed at the reference spread
    q     = setNames(runif(length(conds), -0.2, 0.6), conds),                           # prior position per condition
    lambda = setNames(rep(exp(runif(1, log(0.8), log(3))), length(conds)), conds)       # loss aversion, shared
  )
}

# ---------------------------------------------------------------------------
# Simulate and refit, one dataset at a time (resumes from the checkpoint)
# ---------------------------------------------------------------------------
results <- if (file.exists(cache_file)) readRDS(cache_file) else list()
for (d in seq_len(N_DATASETS)) {
  if (d <= length(results) && !is.null(results[[d]])) next
  cat(sprintf("\n==== dataset %d of %d  %s\n", d, N_DATASETS, format(Sys.time())))
  tp <- true_params[[d]]

  # One simulated choice per real trial, from the closed-form probabilities.
  set.seed(SEED + d)
  pr  <- choice_probabilities(prep, tp$sigma, list(s0 = tp$s0, q = tp$q, lambda = tp$lambda))
  sim <- as.integer(runif(nrow(pr)) < pr$p_left)
  acc <- mean((sim == 1) == (true_diff$Diff > 0), na.rm = TRUE)
  cat(sprintf("  true q = %s  lambda = %.3f  simulated accuracy = %.3f\n",
              paste(sprintf("%.3f", tp$q), collapse = " / "), tp$lambda[[1]], acc))

  # Refit from scratch on the profile likelihood.
  prep_d <- set_choices(prep, sim)
  fit    <- fit_profile(prep_d, q_grid = Q_GRID, joint_maxit = 60, verbose = FALSE)
  cat(sprintf("  recovered q = %s  lambda = %.3f  nll = %.1f\n",
              paste(sprintf("%.3f", fit$q), collapse = " / "), fit$lambda[[1]], fit$nll))

  results[[d]] <- list(
    dataset  = d,
    accuracy = acc,
    nll      = fit$nll,
    sigma = tibble(dataset = d, Study = subjects$Study, Subject = subjects$Subject, WN = subjects$WN,
                   true = unname(tp$sigma), recovered = unname(fit$sigma[subjects$Subject])),
    group = tibble(dataset = d, WN = conds,
                   true_s0 = unname(tp$s0), rec_s0 = unname(fit$s0[conds]),
                   true_q  = unname(tp$q),  rec_q  = unname(fit$q[conds]),
                   true_lambda = unname(tp$lambda[conds]), rec_lambda = unname(fit$lambda[conds]))
  )
  saveRDS(results, cache_file)     # checkpoint
}

# ---------------------------------------------------------------------------
# Summary and figure
# ---------------------------------------------------------------------------
sig_df <- tibble()
grp_df <- tibble()
for (r in results) {
  sig_df <- bind_rows(sig_df, r$sigma)
  grp_df <- bind_rows(grp_df, r$group)
}

cat("\n==== Recovery summary over", length(results), "datasets\n")
for (s in c(1, 2)) {
  x <- sig_df %>% filter(Study == s)
  cat(sprintf("sigma, Study %d: r(log) = %.3f, median |relative error| = %.3f\n", s,
              cor(log(x$true), log(x$recovered)), median(abs(x$recovered - x$true) / x$true)))
}
cat(sprintf("q:  r = %.3f, median |rec - true| = %.3f\n",
            cor(grp_df$true_q, grp_df$rec_q), median(abs(grp_df$rec_q - grp_df$true_q))))
lam_df <- grp_df %>% filter(WN == conds[1])    # lambda is shared, one point per dataset
cat(sprintf("lambda: r(log) = %.3f, median |log(rec / true)| = %.3f\n",
            cor(log(lam_df$true_lambda), log(lam_df$rec_lambda)), median(abs(log(lam_df$rec_lambda / lam_df$true_lambda)))))
print(as.data.frame(grp_df %>% mutate(across(where(is.numeric), \(x) round(x, 3)))))

fig <- plot_parameter_recovery(sig_df, grp_df)
ggsave(file.path(fig_dir, "parameter_recovery_bayesian.pdf"), fig, width = 2 * .figw, height = 2 * .figh)
cat("\nFigure written to", file.path(fig_dir, "parameter_recovery_bayesian.pdf"), "\n")
cat("Done", format(Sys.time()), "\n")
