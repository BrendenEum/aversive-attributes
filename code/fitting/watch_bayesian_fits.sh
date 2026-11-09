#!/bin/zsh
# Live view of the Bayesian fitting queue (run_all_bayesian_fits.sh): the
# running jobs, how many folds and recovery datasets are done, the latest
# result line of the most recently active logs, and memory. Refreshes every
# 30 s; Ctrl-C to stop. From the project root:
#
#   zsh code/fitting/watch_bayesian_fits.sh
cd "$(dirname "$0")/../.." || exit 1
setopt null_glob          # a pattern that matches nothing expands to nothing
L=results/logs
C=results/cache/bayesian_model
while true; do
  clear
  echo "=== Bayesian fitting queue   $(date '+%a %H:%M:%S') ==="
  [ -f $L/bayesian_driver.log ] && cat $L/bayesian_driver.log
  echo
  echo "--- running jobs (elapsed) ---"
  ps -eo etime,command | grep "[e]xec/R .*\(run_bayesian_fits\|recover_bayesian\)" \
    | sed 's| */Library.*--args | |; s| */Library.*recover_bayesian_attributes.R| recovery|'
  echo
  echo "--- done ---"
  for d in $C/loo_*_folds(/); do
    files=($d/fold_*.rds)
    printf "%-14s %2d of 12 folds\n" "$(basename $d | sed 's/loo_//; s/_folds//')" ${#files}
  done
  files=($C/parameter_recovery/dataset_*.rds)
  echo "recovery       ${#files} of 20 datasets"
  echo "full fits      $(for f in $C/fit_*.rds; do basename $f .rds | sed 's/fit_//'; done | tr '\n' ' ')"
  echo "assembled      $(for f in $C/loo_*.rds; do basename $f .rds | sed 's/loo_//'; done | tr '\n' ' ')"
  echo
  echo "--- latest result line of the most recently active logs ---"
  for f in $L/(bayesian_fits_*|recovery_*).log(om[1,8]); do     # om = newest first
    line=$(grep -h "^fold\|recovered q\|: nll = " $f | tail -1 | cut -c1-140)
    [ -n "$line" ] && echo "$(basename $f .log): $line"
  done
  echo
  top -l 1 -n 0 2>/dev/null | grep "^PhysMem"
  sleep 30
done
