# Fit every variant of the Bayesian attribute-level model to the real data and
# run the trial-level cross-validation, outside of RStudio. From the project root:
#
#   caffeinate -i nohup Rscript code/fitting/run_bayesian_fits.R > results/logs/bayesian_fits.log 2>&1 &
#
# Arguments: variant names (default: all of them); `fit` = only the full fit;
# `loo` = only the cross-validation; `folds=a-b` = only folds a to b of it
# (`folds=k` for one fold). Every fold is a full refit and is cached in its own
# file (loo_<variant>_folds/fold_k.rds), so a variant's 12 folds can be split
# over processes; whichever run finds all 12 files assembles loo_<variant>.rds.
# With the exact choice probabilities a main-model fit takes about an hour, so
# for example the full fit plus one process per fold:
#
#   nohup Rscript code/fitting/run_bayesian_fits.R cond_lam fit > results/logs/bayesian_fits_cond_lam_fit.log 2>&1 &
#   for f in $(seq 1 12); do
#     nohup Rscript code/fitting/run_bayesian_fits.R cond_lam loo folds=$f > results/logs/bayesian_fits_cond_lam_fold$f.log 2>&1 &
#   done
#
# code/fitting/run_all_bayesian_fits.sh runs every fit, fold and recovery
# dataset as one queue.
#
# Caches go to results/cache/bayesian_model/ (fit_<variant>.rds and
# loo_<variant>.rds); notebooks/bayesian_model.Rmd loads them. Fits and folds
# whose caches exist are skipped, so a run can be repeated after an interruption.

suppressPackageStartupMessages(library(tidyverse))
source(here::here("code/utilities/getAllUtilities.R"))          # paths, build_pas(), build_boxlong()
source(here::here("code/simulation/simulate_bayesian_attributes.R"))
source(here::here("code/fitting/fit_bayesian_attributes.R"))

cache_dir <- file.path(paths$root, paths$results$cache, "bayesian_model")

pas_raw <- read_csv(file.path(paths$root, paths$data$processed, "main_data.csv"), show_col_types = FALSE)
pas     <- build_pas(pas_raw)
boxlong <- build_boxlong(pas)

# Command line: what to run, which folds, which variants.
args  <- commandArgs(trailingOnly = TRUE)
what  <- if ("fit" %in% args) "fit" else if ("loo" %in% args) "loo" else "both"
fspec <- grep("^folds=", args, value = TRUE)
folds <- NULL
if (length(fspec) > 0) {
  r     <- as.integer(strsplit(sub("^folds=", "", fspec[1]), "-")[[1]])
  folds <- seq(r[1], r[length(r)])
}
to_run <- setdiff(args, c("fit", "loo", fspec))
if (length(to_run) == 0) to_run <- names(bayesian_variants())
stopifnot(all(to_run %in% names(bayesian_variants())))

for (name in to_run) {
  run <- run_variant(name, boxlong, cache_dir, what = what, folds = folds)
  if (!is.null(run$fit)) {
    cat(sprintf("\n%s: nll = %.1f  s0 = %s  q = %s  lambda = %s  exp(mu) = %s  tau = %s\n", name, run$fit$nll,
                paste(sprintf("%.3f", run$fit$s0), collapse = " / "),
                paste(sprintf("%.3f", run$fit$q),  collapse = " / "),
                paste(sprintf("%.3f", run$fit$lambda), collapse = " / "),
                paste(sprintf("%.2f", exp(run$fit$mu)), collapse = " / "),
                paste(sprintf("%.2f", run$fit$tau), collapse = " / ")))
  }
  if (!is.null(run$loo)) {
    cat(sprintf("%s: cross-validation assembled; folds not converged: %d (inner: %d)\n", name,
                sum(run$loo$fold_params$convergence != 0), sum(run$loo$fold_params$convergence_inner != 0)))
  }
}
cat("\nDone", format(Sys.time()), "\n")
