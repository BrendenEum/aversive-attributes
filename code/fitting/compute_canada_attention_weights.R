# Fit the main text's attention-weight regression (brms, hierarchical, the
# exact specification and sampler settings of the data fits in
# Main_Analysis.Rmd) to the choices simulated from each Bayesian model
# variant, on a Compute Canada (Digital Research Alliance) cluster.
#
# One hundred fits: five variants (per-condition prior with and without loss
# aversion, per-type prior, zero prior, flat prior) x two range conditions x
# ten replicates of the simulated choices. A replicate is one simulated
# choice per trial, drawn from the variant's held-out probabilities; the
# replicates are fit as separate datasets (each the size of the data, so the
# hierarchical model sees the data's design) and the notebooks pool their
# posteriors. Each fit is one SLURM array task (four chains, one core each,
# plus within-chain threading over the task's remaining cores). The fits give
# the model side of the attention-weight overlays and of the OOS-MAE column
# of tables/bayesian_model_fits.tex on the same estimator as the data side.
#
# Workflow (the exact commands are in compute_canada_attention_weights.sh):
#
#   1. Locally, from the project root. Draws the ten replicates of simulated
#      choices from each variant's held-out (12-fold) probabilities (seed 1;
#      replicate k is the k-th of the notebook's ten) and writes the
#      regression's design matrix per variant, range and replicate to
#      results/compute_canada/inputs/sim_<variant>_<range>_rep<k>.rds:
#
#        Rscript code/fitting/compute_canada_attention_weights.R prepare
#
#   2. Copy this script, the .sh and the inputs into one folder on the
#      cluster (the inputs are not in git: *.rds is ignored and they carry the
#      participants' dwell times). On the login node, once,
#
#        Rscript compute_canada_attention_weights.R setup
#
#      installs brms, posterior, dplyr, cmdstanr and CmdStan into the user
#      library if they are missing (compute nodes have no internet). Then the
#      array job runs, per task,
#
#        Rscript compute_canada_attention_weights.R fit <1-100>
#
#      which leaves the full brmsfit in outputs/ (about 1.4 GB, stays on the
#      cluster) and a compact file next to it,
#      attention_weights_utility_sim_<variant>_<range>_rep<k>_draws.rds: the
#      population-level draws, their summaries, and the sampler diagnostics.
#
#   3. Copy the compact files back into results/compute_canada/outputs/ and
#      check them; they are moved into notebooks/brm_fits/, where
#      bayesian_model.Rmd reads them:
#
#        Rscript code/fitting/compute_canada_attention_weights.R collect
#
# Other modes: `list` prints the job table (the array index of every
# variant x range x replicate); `fit <i> --smoke` is a two-minute test of the fit mode
# (eight subjects, one short chain) that writes under smoke/ and touches
# nothing else.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "list" else args[1]

# Two layouts, told apart by the script's own location. In the project
# (code/fitting/) the inputs and outputs live under results/compute_canada/.
# On the cluster the script sits in one flat folder next to the input files,
# and writes outputs/ and smoke/ inside that folder.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir  <- if (length(script_path) == 1) normalizePath(dirname(script_path)) else getwd()
in_project  <- file.exists(file.path(script_dir, "..", "..", "notebooks", "bayesian_model.Rmd"))
root        <- if (in_project) normalizePath(file.path(script_dir, "..", "..")) else script_dir

in_dir    <- if (in_project) file.path(root, "results", "compute_canada", "inputs")  else root
out_dir   <- if (in_project) file.path(root, "results", "compute_canada", "outputs") else file.path(root, "outputs")
smoke_dir <- if (in_project) file.path(root, "results", "compute_canada", "smoke")   else file.path(root, "smoke")
fits_dir  <- file.path(root, "notebooks", "brm_fits")

# The jobs: replicate 1 of every variant x range first (the main model's two
# ranges, then the four supplementary variants), then replicate 2, and so on.
# The SLURM array index is the row number.
REPS <- 10
JOBS <- data.frame(
  rep     = rep(seq_len(REPS), each = 10),
  variant = rep(rep(c("cond_lam", "cond_lam1", "type_lam", "zero", "flat"), each = 2), times = REPS),
  range   = rep(c("Narrow", "Wide"), times = 5 * REPS),
  stringsAsFactors = FALSE
)
JOBS$label <- c(cond_lam  = "Per-condition prior",
                cond_lam1 = "Per-condition prior, no loss aversion",
                type_lam  = "Per-type prior",
                zero      = "Zero prior",
                flat      = "Flat prior")[JOBS$variant]
job_stem   <- function(v, rng, k) paste0(v, "_", tolower(rng), "_rep", k)
input_file <- function(v, rng, k) file.path(in_dir, paste0("sim_", job_stem(v, rng, k), ".rds"))
draws_name <- function(v, rng, k) paste0("attention_weights_utility_sim_", job_stem(v, rng, k), "_draws.rds")

# Sampler settings: the same as the data fits (Main_Analysis.Rmd, chunk
# flex-fit-by-range) so that both sides of the overlay are the same estimator.
SETTINGS <- list(chains = 4, iter = 4000, seed = 1234, adapt_delta = 0.95, max_treedepth = 12,
                 prior_b = "normal(0, 2.5)", prior_intercept = "normal(0, 2)", backend = "cmdstanr")

# ---------------------------------------------------------------------------
# list
# ---------------------------------------------------------------------------
if (mode == "list") {
  cat("Array index, replicate, variant, range:\n")
  print(cbind(index = seq_len(nrow(JOBS)), JOBS), row.names = FALSE)
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# setup (cluster): the R packages and CmdStan, into the user library, once.
# ---------------------------------------------------------------------------
if (mode == "setup") {
  lib <- path.expand(Sys.getenv("R_LIBS_USER"))        # R's default personal library for this R version
  if (!nzchar(lib)) stop("R_LIBS_USER is not set.")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2")))
  Sys.setenv(MAKEFLAGS = paste0("-j", cores))
  # An interrupted install leaves 00LOCK-<package> folders that block the next
  # one; an interrupted CmdStan build leaves a folder without the example
  # binary that install_cmdstan() compiles last. Clear both before checking.
  locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE, include.dirs = TRUE)
  if (length(locks)) { cat("Removing stale install locks:", basename(locks), "\n"); unlink(locks, recursive = TRUE) }
  need    <- c("brms", "posterior", "dplyr", "cmdstanr")
  missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    cat("Installing", paste(missing, collapse = ", "), "into", lib, "\n")
    install.packages(missing, lib = lib, repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
    still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still)) stop("Failed to install: ", paste(still, collapse = ", "))
  }
  cmdstan_complete <- tryCatch({
    d <- cmdstanr::cmdstan_path()
    file.exists(file.path(d, "bin", "stanc")) && length(list.files(file.path(d, "examples", "bernoulli"), pattern = "^bernoulli$|^bernoulli\\.exe$")) > 0
  }, error = function(e) FALSE)
  if (!cmdstan_complete) {
    partial <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
    if (!is.null(partial)) { cat("Removing incomplete CmdStan build:", partial, "\n"); unlink(partial, recursive = TRUE) }
    cmdstanr::install_cmdstan(cores = cores)
  }
  cat(sprintf("Setup ok: R %s, brms %s, CmdStan %s, library %s\n", getRversion(),
              packageVersion("brms"), cmdstanr::cmdstan_version(), lib))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# prepare (local): ten replicates of one simulated choice per trial from each
# variant's held-out probabilities, through the notebook's own design builder.
# ---------------------------------------------------------------------------
if (mode == "prepare") {
  rmd    <- file.path(root, "notebooks", "bayesian_model.Rmd")
  lines  <- readLines(rmd); starts <- grep("^```\\{r [^,}]+", lines); ends <- grep("^```$", lines)
  chunks <- list()
  for (i in starts) {
    label <- sub("^```\\{r ([^,}]+).*$", "\\1", lines[i])
    chunks[[label]] <- lines[(i + 1):(min(ends[ends > i]) - 1)]
  }
  run_chunk <- function(label) {
    cat("\n==== chunk:", label, "====", format(Sys.time()), "\n")
    eval(parse(text = chunks[[label]]), envir = globalenv())
  }
  run_chunk("libraries"); run_chunk("options")
  RUN_FITS <- FALSE                                   # the variants' caches only
  for (ch in c("plot-helpers", "read-data", "flex-design", "model-source", "model-helpers", "fits")) run_chunk(ch)

  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- data.frame()
  for (v in unique(JOBS$variant)) {
    if (is.null(runs[[v]])) stop("No cached fit for variant ", v, " in results/cache/bayesian_model/.")
    # The notebook's ten replicates (seed 1): replicate k is the k-th draw.
    tr       <- runs[[v]]$loo$pred %>% transmute(Subject, Trial, WN, Study, p_left = p_left_heldout)
    sims_all <- simulate_choices_from_prob(tr, reps = REPS, seed = 1)
    for (k in seq_len(REPS)) {
      sims <- sims_all %>% filter(rep == k) %>%
        transmute(Subject, Trial, WN, Study, TrialKey = Trial, ChooseL = ChooseL_sim)
      bl_sim <- boxlong %>% select(-ChooseL) %>%
        inner_join(sims %>% select(Subject, Trial, TrialKey, ChooseL), by = c("Subject", "Trial"))
      for (rng in c("Narrow", "Wide")) {
        dat <- build_utility_dat(rng, bl_sim, sims) %>% rename(Trial = TrialKey)
        stopifnot(!anyNA(dat), all(c("Subject", "Trial", "ChooseL", "S2") %in% names(dat)))
        saveRDS(dat, input_file(v, rng, k))
        manifest <- rbind(manifest, data.frame(
          variant = v, range = rng, rep = k, file = basename(input_file(v, rng, k)), n_trials = nrow(dat),
          n_subjects = n_distinct(dat$Subject), mean_choose_left = round(mean(dat$ChooseL), 4),
          n_utility_terms = sum(grepl("^U_", names(dat))), n_dwell_terms = sum(grepl("^[PNB]_(neg|pos)[0-9]+$", names(dat))),
          sim_seed = 1, prepared = format(Sys.time())))
        cat(sprintf("%-10s %-6s rep %2d: %6d trials, %3d subjects, %2d utilities, %2d cells -> %s\n",
                    v, rng, k, nrow(dat), n_distinct(dat$Subject), sum(grepl("^U_", names(dat))),
                    sum(grepl("^[PNB]_(neg|pos)[0-9]+$", names(dat))), basename(input_file(v, rng, k))))
      }
    }
  }
  write.csv(manifest, file.path(in_dir, "manifest.csv"), row.names = FALSE)
  cat("\nInputs written to", in_dir, "\n")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# fit (cluster): one job = one variant x range x replicate.
# ---------------------------------------------------------------------------
if (mode == "fit") {
  if (length(args) < 2) stop("Usage: fit <array index 1-100 | variant range replicate> [--smoke]")
  smoke <- "--smoke" %in% args
  if (grepl("^[0-9]+$", args[2])) {
    i <- as.integer(args[2]); v <- JOBS$variant[i]; rng <- JOBS$range[i]; k <- JOBS$rep[i]
  } else {
    v <- args[2]; rng <- args[3]; k <- as.integer(args[4])
  }
  if (is.na(v) || !(v %in% JOBS$variant) || !(rng %in% c("Narrow", "Wide")) || is.na(k) || !(k %in% seq_len(REPS))) {
    stop("Unknown job: ", paste(args[-1], collapse = " "))
  }

  suppressPackageStartupMessages({ library(brms); library(posterior); library(dplyr) })
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
  options(brms.backend = SETTINGS$backend, mc.cores = cores)

  dat <- readRDS(input_file(v, rng, k))
  util_cols  <- grep("^U_", names(dat), value = TRUE)
  dwell_cols <- grep("^[PNB]_(neg|pos)[0-9]+$", names(dat), value = TRUE)
  # build_flex_formula() of Main_Analysis.Rmd: population utility for every
  # value and population weight for every cell, subject-level random slopes on
  # every utility (uncorrelated) and on every cell weight (correlated).
  formula <- as.formula(paste(
    "ChooseL ~ 1 + S2 +", paste(util_cols, collapse = " + "), "+", paste(dwell_cols, collapse = " + "),
    "+ (0 +", paste(util_cols,  collapse = " + "), "|| Subject)",
    "+ (0 +", paste(dwell_cols, collapse = " + "), "| Subject)"
  ))
  priors <- c(prior(normal(0, 2.5), class = "b"), prior(normal(0, 2), class = "Intercept"))

  chains <- SETTINGS$chains; iter <- SETTINGS$iter
  # Cores beyond one per chain go to within-chain threading (reduce_sum): the
  # same posterior, a faster fit. Four cores = four plain chains.
  threads <- max(1L, cores %/% chains)
  write_dir <- out_dir
  if (smoke) {
    keep <- head(unique(dat$Subject), 8)
    dat  <- dat %>% filter(Subject %in% keep)
    chains <- 1; iter <- 300; write_dir <- smoke_dir
  }
  dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
  fit_file   <- file.path(write_dir, paste0("fit_", job_stem(v, rng, k)))        # brm adds .rds
  draws_file <- file.path(write_dir, draws_name(v, rng, k))

  cat(sprintf("Fitting %s (%s), %s, replicate %d: %d trials, %d subjects, %d cores, %d chains x %d iterations, %d thread(s) per chain%s\n",
              v, JOBS$label[JOBS$variant == v][1], rng, k, nrow(dat), n_distinct(dat$Subject), cores, chains, iter, threads,
              if (smoke) " [SMOKE TEST]" else ""))
  t0  <- Sys.time()
  fit <- brm(formula, data = dat, family = bernoulli(), prior = priors,
             chains = chains, cores = min(chains, cores), iter = iter, seed = SETTINGS$seed,
             control = list(adapt_delta = SETTINGS$adapt_delta, max_treedepth = SETTINGS$max_treedepth),
             threads = threading(threads), file = fit_file, file_refit = "never")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))

  # The compact file: population-level draws (the intercepts, the utilities,
  # the attention weights), their summaries and the sampler diagnostics.
  draws <- as_draws_df(fit, variable = "^b_", regex = TRUE)
  summ  <- summarise_draws(draws, "mean", "sd", ~quantile2(.x, probs = c(0.025, 0.975)), "rhat", "ess_bulk", "ess_tail")
  np    <- nuts_params(fit)
  out <- list(
    variant = v, label = JOBS$label[JOBS$variant == v][1], range = rng, rep = k, smoke = smoke,
    draws = draws, summary = summ,
    n_divergent     = sum(np$Value[np$Parameter == "divergent__"]),
    n_max_treedepth = sum(np$Value[np$Parameter == "treedepth__"] >= SETTINGS$max_treedepth),
    n_trials = nrow(dat), n_subjects = n_distinct(dat$Subject),
    formula = deparse1(formula), settings = modifyList(SETTINGS, list(chains = chains, iter = iter, threads_per_chain = threads)),
    versions = list(R = R.version.string, brms = as.character(packageVersion("brms")),
                    cmdstan = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA)),
    elapsed_hours = elapsed, finished = Sys.time()
  )
  saveRDS(out, draws_file)
  cat(sprintf("Done in %.2f h: %d divergent transitions, %d max-treedepth hits, max Rhat %.3f, min bulk ESS %.0f\n%s\n",
              elapsed, out$n_divergent, out$n_max_treedepth, max(summ$rhat, na.rm = TRUE), min(summ$ess_bulk, na.rm = TRUE), draws_file))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# collect (local): check the compact files that came back and move them to
# notebooks/brm_fits/, where bayesian_model.Rmd reads them.
# ---------------------------------------------------------------------------
if (mode == "collect") {
  dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- data.frame()
  for (i in seq_len(nrow(JOBS))) {
    v <- JOBS$variant[i]; rng <- JOBS$range[i]; k <- JOBS$rep[i]
    src <- file.path(out_dir, draws_name(v, rng, k)); dst <- file.path(fits_dir, draws_name(v, rng, k))
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
    if (!file.exists(dst)) { rows <- rbind(rows, data.frame(variant = v, range = rng, rep = k, status = "missing")); next }
    x <- readRDS(dst)
    w <- x$summary[grepl("^b_[PNB]_(neg|pos)[0-9]+$", x$summary$variable), ]
    rows <- rbind(rows, data.frame(
      variant = v, range = rng, rep = k, status = if (isTRUE(x$smoke)) "SMOKE" else "ok",
      hours = round(x$elapsed_hours, 2), divergent = x$n_divergent, treedepth = x$n_max_treedepth,
      max_rhat = round(max(x$summary$rhat, na.rm = TRUE), 3), min_ess_bulk = round(min(x$summary$ess_bulk, na.rm = TRUE)),
      weights_min = round(min(w$mean), 3), weights_max = round(max(w$mean), 3)))
  }
  print(rows, row.names = FALSE)
  # One line per variant x range: how many replicates are back, and the worst diagnostics over them.
  ok <- rows[rows$status == "ok", ]
  if (nrow(ok)) {
    cat("\nPer variant x range (replicates back, hours per fit, worst diagnostics):\n")
    agg <- do.call(rbind, lapply(split(ok, list(ok$variant, ok$range), drop = TRUE), function(d) data.frame(
      variant = d$variant[1], range = d$range[1], reps = nrow(d), hours_mean = round(mean(d$hours), 2), hours_max = round(max(d$hours), 2),
      divergent = sum(d$divergent), treedepth = sum(d$treedepth), max_rhat = max(d$max_rhat), min_ess_bulk = min(d$min_ess_bulk))))
    print(agg[order(match(agg$variant, unique(JOBS$variant)), agg$range), ], row.names = FALSE)
  }
  quit(save = "no")
}

stop("Unknown mode '", mode, "'. Use list, prepare, setup, fit <i> [--smoke], or collect.")
