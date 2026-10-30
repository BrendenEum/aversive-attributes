# Fit every variant of the Bayesian attribute-level model to the real data and
# run the trial-level leave-one-out, outside of RStudio. From the project root:
#
#   caffeinate -i nohup Rscript code/fitting/run_bayesian_fits.R > results/logs/bayesian_fits.log 2>&1 &
#
# Variant names on the command line restrict the run to those variants, so the
# work can be split over several processes (every fold is a full refit, about
# 2.5 minutes each for the main model and about a minute with lambda or q
# fixed, 60 folds per variant; all variants took about 3.5 hours in four
# processes on 2026-10-30, with the main model alone in one of them):
#
#   nohup Rscript code/fitting/run_bayesian_fits.R cond_lam > results/logs/bayesian_fits_cond_lam.log 2>&1 &
#
# Caches go to results/cache/bayesian_model/ (fit_<variant>.rds and
# loo_<variant>.rds); notebooks/bayesian_model.Rmd loads them. Variants whose
# caches exist are skipped, so the script can be re-run after an interruption.

suppressPackageStartupMessages(library(tidyverse))
source(here::here("code/utilities/getAllUtilities.R"))          # paths, build_pas(), build_boxlong()
source(here::here("code/simulation/simulate_bayesian_attributes.R"))
source(here::here("code/fitting/fit_bayesian_attributes.R"))

cache_dir <- file.path(paths$root, paths$results$cache, "bayesian_model")

pas_raw <- read_csv(file.path(paths$root, paths$data$processed, "main_data.csv"), show_col_types = FALSE)
pas     <- build_pas(pas_raw)
boxlong <- build_boxlong(pas)

# Which variants to run: those named on the command line, or all of them.
to_run <- commandArgs(trailingOnly = TRUE)
if (length(to_run) == 0) to_run <- names(bayesian_variants())
stopifnot(all(to_run %in% names(bayesian_variants())))

for (name in to_run) {
  run <- run_variant(name, boxlong, cache_dir)
  cat(sprintf("\n%s: nll = %.1f  s0 = %s  q = %s  lambda = %s  folds not converged: %d\n", name, run$fit$nll,
              paste(sprintf("%.3f", run$fit$s0), collapse = " / "),
              paste(sprintf("%.3f", run$fit$q),  collapse = " / "),
              paste(sprintf("%.3f", run$fit$lambda), collapse = " / "),
              sum(run$loo$fold_params$convergence != 0)))
}
cat("\nDone", format(Sys.time()), "\n")
