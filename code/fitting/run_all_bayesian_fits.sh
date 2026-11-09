#!/bin/zsh
# Run every Bayesian-model fit, cross-validation fold and parameter-recovery
# dataset as a queue of jobs, seven at a time, then assemble the results.
# From the project root:
#
#   caffeinate -i nohup zsh code/fitting/run_all_bayesian_fits.sh > results/logs/bayesian_driver.log 2>&1 &
#
# Every job caches its own output (see run_bayesian_fits.R and
# recover_bayesian_attributes.R), so the script can be re-run after an
# interruption and only the missing pieces are computed. Every fold and every
# recovery dataset is a full fit from scratch. With the exact choice
# probabilities one takes about half an hour alone and longer when seven jobs
# share the machine (79 min measured on 2026-10-30, on an 8-core M2 with idle
# apps holding most of the memory). Jobs are one fold or one dataset each,
# so the seven slots stay busy to the end; the whole queue is roughly 100
# process-hours.
cd "$(dirname "$0")/../.." || exit 1
R=code/fitting/run_bayesian_fits.R
REC=code/parameter_recovery/recover_bayesian_attributes.R
L=results/logs
mkdir -p "$L"

jobs=()
# Main model first: the full fit (a no-op when cached), then its twelve folds.
jobs+=("Rscript $R cond_lam fit > $L/bayesian_fits_cond_lam_fit.log 2>&1")
for f in $(seq 1 12); do
  jobs+=("Rscript $R cond_lam loo folds=$f > $L/bayesian_fits_cond_lam_fold$f.log 2>&1")
done
# Parameter recovery, one dataset per job.
for d in $(seq 1 20); do
  jobs+=("Rscript $REC datasets=$d-$d > $L/recovery_$d.log 2>&1")
done
# The other variants that fit lambda: the full fit, then twelve folds each.
for v in type_lam cond_lamcond zero flat; do
  jobs+=("Rscript $R $v fit > $L/bayesian_fits_${v}_fit.log 2>&1")
  for f in $(seq 1 12); do
    jobs+=("Rscript $R $v loo folds=$f > $L/bayesian_fits_${v}_fold$f.log 2>&1")
  done
done
# The lambda = 1 variants are cheap (the probit is exact there): one job each.
jobs+=("Rscript $R cond_lam1 > $L/bayesian_fits_cond_lam1.log 2>&1")
jobs+=("Rscript $R type_lam1 > $L/bayesian_fits_type_lam1.log 2>&1")

echo "queue of ${#jobs[@]} jobs, seven at a time, started $(date)"
printf '%s\n' "${jobs[@]}" | xargs -P 7 -I CMD sh -c 'CMD'
echo "all jobs finished $(date)"

# Assemble whatever is not yet assembled and print the summaries.
Rscript $R > $L/bayesian_fits_assemble.log 2>&1
Rscript $REC > $L/recovery_summary.log 2>&1
echo "Done $(date)"
