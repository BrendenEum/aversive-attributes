# Parameter recovery for the main Bayesian attribute-level model (variant
# cond_lam: one prior per range condition, shared loss aversion, per-subject
# noise a log-normal random effect with its mean and SD fitted per range
# condition; see code/fitting/fit_bayesian_attributes.R).
#
# From the project root:
#
#   caffeinate -i nohup Rscript code/parameter_recovery/recover_bayesian_attributes.R > results/logs/recovery.log 2>&1 &
#
# For each of N_DATASETS datasets (20; each takes about an hour with the exact
# choice probabilities): draw true parameters (the subjects' sigmas from the
# true population distribution of their condition), simulate one choice per
# real trial (real stimuli, real dwell times), refit from scratch with the
# profile-likelihood procedure in code/fitting/fit_bayesian_attributes.R, and
# store true and recovered values (a subject's recovered sigma is the posterior
# mean). Trials with an unhovered box are dropped, as in the real fits.
# Every dataset is cached in its own file
# (results/cache/bayesian_model/parameter_recovery/dataset_dd.rds), so the
# datasets can be split over processes with `datasets=a-b` on the command line:
#
#   for d in 1-4 5-8 9-12 13-16 17-20; do
#     nohup Rscript code/parameter_recovery/recover_bayesian_attributes.R datasets=$d > results/logs/recovery_$d.log 2>&1 &
#   done
#
# Whichever run finds all N_DATASETS files writes the combined
# parameter_recovery.rds (which the notebook reads), prints the summary and
# draws the figure (plot_parameter_recovery() in code/utilities/) to
# paper/figures/parameter_recovery_bayesian.pdf. Re-running with a complete
# cache only redraws the figure.

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

variant  <- bayesian_variants()$cond_lam      # per-condition prior, shared loss aversion
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
  tp <- list(
    mu     = setNames(runif(length(conds), log(4), log(12)), conds),                    # population mean of log sigma per condition (fitted: exp(mu) ~ 6 Narrow, 9 Wide points per second)
    tau    = setNames(runif(length(conds), 0.3, 0.8), conds),                           # population SD of log sigma per condition (fitted ~0.5)
    s0     = setNames(rep(1, length(conds)), conds),                                    # prior SD fixed at the reference spread
    q      = setNames(runif(length(conds), -0.2, 0.6), conds),                          # prior position per condition
    lambda = setNames(rep(exp(runif(1, log(0.8), log(3))), length(conds)), conds)       # loss aversion, shared
  )
  # Noise per subject: a draw from their condition's population distribution.
  wn <- as.character(subjects$WN)
  tp$sigma <- setNames(exp(rnorm(nrow(subjects), tp$mu[wn], tp$tau[wn])), subjects$Subject)
  true_params[[d]] <- tp
}

# ---------------------------------------------------------------------------
# Simulate and refit, one dataset at a time (each cached in its own file)
# ---------------------------------------------------------------------------
# Which datasets this process does: datasets=a-b on the command line, or all.
args  <- commandArgs(trailingOnly = TRUE)
dspec <- grep("^datasets=", args, value = TRUE)
todo  <- seq_len(N_DATASETS)
if (length(dspec) > 0) {
  r    <- as.integer(strsplit(sub("^datasets=", "", dspec[1]), "-")[[1]])
  todo <- seq(r[1], r[length(r)])
}
set_dir  <- file.path(cache_dir, "parameter_recovery")
set_file <- function(d) file.path(set_dir, sprintf("dataset_%02d.rds", d))
dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)

for (d in todo) {
  if (file.exists(set_file(d))) next
  cat(sprintf("\n==== dataset %d of %d  %s\n", d, N_DATASETS, format(Sys.time())))
  tp <- true_params[[d]]

  # One simulated choice per real trial, from the model's choice probabilities.
  set.seed(SEED + d)
  pr  <- choice_probabilities(prep, tp$sigma, list(s0 = tp$s0, q = tp$q, lambda = tp$lambda))
  sim <- as.integer(runif(nrow(pr)) < pr$p_left)
  acc <- mean((sim == 1) == (true_diff$Diff > 0), na.rm = TRUE)
  cat(sprintf("  true q = %s  lambda = %.3f  exp(mu) = %s  tau = %s  simulated accuracy = %.3f\n",
              paste(sprintf("%.3f", tp$q), collapse = " / "), tp$lambda[[1]],
              paste(sprintf("%.2f", exp(tp$mu)), collapse = " / "), paste(sprintf("%.2f", tp$tau), collapse = " / "), acc))

  # Refit from scratch with the same three stages as the real data (only the
  # stage-1 grid is coarser).
  prep_d <- set_choices(prep, sim)
  fit    <- fit_profile(prep_d, q_grid = Q_GRID, verbose = FALSE)
  cat(sprintf("  recovered q = %s  lambda = %.3f  exp(mu) = %s  tau = %s  nll = %.1f  convergence = %s / %d\n",
              paste(sprintf("%.3f", fit$q), collapse = " / "), fit$lambda[[1]],
              paste(sprintf("%.2f", exp(fit$mu)), collapse = " / "), paste(sprintf("%.2f", fit$tau), collapse = " / "),
              fit$nll, paste(fit$convergence, collapse = " "), fit$convergence_inner))

  saveRDS(list(
    dataset  = d,
    accuracy = acc,
    nll      = fit$nll,
    sigma = tibble(dataset = d, Study = subjects$Study, Subject = subjects$Subject, WN = subjects$WN,
                   true = unname(tp$sigma), recovered = unname(fit$sigma[subjects$Subject])),
    group = tibble(dataset = d, WN = conds,
                   true_s0 = unname(tp$s0), rec_s0 = unname(fit$s0[conds]),
                   true_q  = unname(tp$q),  rec_q  = unname(fit$q[conds]),
                   true_lambda = unname(tp$lambda[conds]), rec_lambda = unname(fit$lambda[conds]),
                   true_mu  = unname(tp$mu[conds]),  rec_mu  = unname(fit$mu[conds]),
                   true_tau = unname(tp$tau[conds]), rec_tau = unname(fit$tau[conds]))
  ), set_file(d))
}

# The summary and figure need every dataset; another process may still be on some.
done <- file.exists(sapply(seq_len(N_DATASETS), set_file))
if (!all(done)) {
  cat(sprintf("\n%d of %d datasets done; the summary and figure are made once all exist. %s\n", sum(done), N_DATASETS, format(Sys.time())))
  quit(save = "no")
}
results <- lapply(seq_len(N_DATASETS), function(d) readRDS(set_file(d)))
saveRDS(results, cache_file)

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
cat(sprintf("mu:  r = %.3f, median |rec - true| = %.3f\n",
            cor(grp_df$true_mu, grp_df$rec_mu), median(abs(grp_df$rec_mu - grp_df$true_mu))))
cat(sprintf("tau: r = %.3f, median |rec - true| = %.3f\n",
            cor(grp_df$true_tau, grp_df$rec_tau), median(abs(grp_df$rec_tau - grp_df$true_tau))))
print(as.data.frame(grp_df %>% mutate(across(where(is.numeric), \(x) round(x, 3)))))

fig <- plot_parameter_recovery(sig_df, grp_df)
ggsave(file.path(fig_dir, "parameter_recovery_bayesian.pdf"), fig, width = 2 * .figw, height = 3 * .figh)
cat("\nFigure written to", file.path(fig_dir, "parameter_recovery_bayesian.pdf"), "\n")
cat("Done", format(Sys.time()), "\n")
