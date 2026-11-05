#!/bin/bash
# SLURM array job: the ten brms attention-weight fits on the Bayesian model
# variants' simulated choices (see compute_canada_attention_weights.R), one
# variant x range per task. Four chains, one core each, plus within-chain
# threading (reduce_sum) with the cores beyond four: 8 cores = 2 threads per
# chain, about 1.3-1.7x faster; the same posterior either way.
#
# Written for Fir (fir.alliancecan.ca): 192-core AMD EPYC 9655 nodes with
# 768 GB, i.e. 4 GB per core. Each task asks for 16 cores (four chains x four
# threads) and 32 GB. Set ACCOUNT below to your allocation (the def-... folder
# in ~/projects on the cluster) before submitting.
#
# One-time setup, on the Fir LOGIN node (compute nodes have no internet, so
# packages cannot be installed inside a job). Installing brms compiles rstan
# and building CmdStan takes a while: expect 30-60 minutes in all.
#
#   module load StdEnv/2023 gcc/12.3 r/4.4.0      # `module spider r` lists the R versions
#   mkdir -p ~/.local/R/$EBVERSIONR && export R_LIBS=~/.local/R/$EBVERSIONR
#   export MAKEFLAGS=-j4
#   Rscript -e 'install.packages(c("brms", "posterior", "dplyr"), repos = "https://cloud.r-project.org")'
#   Rscript -e 'install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))'
#   Rscript -e 'cmdstanr::install_cmdstan(cores = 4)'
#   Rscript -e 'library(brms); cat("brms", as.character(packageVersion("brms")), "cmdstan", as.character(cmdstanr::cmdstan_version()), "\n")'
#
# Getting the files there, from the project root on your laptop (the inputs
# are not in git; -R keeps the paths, which the R script relies on):
#
#   rsync -avR code/fitting/compute_canada_attention_weights.R code/fitting/compute_canada_attention_weights.sh \
#         results/compute_canada/inputs/ USER@fir.alliancecan.ca:~/scratch/aversive-attributes/
#
# Submitting, from ~/scratch/aversive-attributes on Fir. First a 20-minute
# smoke test of one task (checks the modules, the library and the compiler
# inside a job), then the real array:
#
#   mkdir -p results/logs
#   sbatch --array=1 --time=0:20:00 --cpus-per-task=4 --mem=8G --export=ALL,SMOKE=1 code/fitting/compute_canada_attention_weights.sh
#   squeue -u $USER                            # wait for it to finish, then read its log in results/logs/
#   sbatch code/fitting/compute_canada_attention_weights.sh
#
# Bringing the compact files back, from the project root on your laptop (the
# full 1.4 GB fits stay on the cluster):
#
#   rsync -av 'USER@fir.alliancecan.ca:~/scratch/aversive-attributes/results/compute_canada/outputs/*_draws.rds' results/compute_canada/outputs/
#   rsync -av 'USER@fir.alliancecan.ca:~/scratch/aversive-attributes/results/logs/compute_canada_*.log' results/logs/
#   Rscript code/fitting/compute_canada_attention_weights.R collect
#
# Time: the same fit on the data took 6.6 hours on four M2 cores with one
# thread per chain. Four threads per chain should bring a Fir task to a few
# hours; the 24-hour limit is a safety margin (a task that hits the limit
# loses its fit, since brms writes the file at the end).

#SBATCH --job-name=attn_sim
#SBATCH --account=def-CHANGEME      # your allocation: the def-... folder in ~/projects
#SBATCH --array=1-10
#SBATCH --cpus-per-task=16          # 4 chains x 4 threads per chain (4 cores = plain chains)
#SBATCH --mem=32G                   # within the 4 GB per core of a Fir node
#SBATCH --time=1-00:00
#SBATCH --output=results/logs/compute_canada_%A_%a.log
##SBATCH --mail-type=END,FAIL       # uncomment both lines for an email when a task ends
##SBATCH --mail-user=you@example.com

module load StdEnv/2023 gcc/12.3 r/4.4.0      # the same modules as the setup
export R_LIBS=~/.local/R/$EBVERSIONR          # the personal library the setup installed into

export TMPDIR="$SLURM_TMPDIR"          # CmdStan's per-chain CSV output goes to node-local disk
export OMP_NUM_THREADS=1               # Stan threads use TBB, not OpenMP

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "Task $SLURM_ARRAY_TASK_ID on $(hostname), $(date)"
Rscript code/fitting/compute_canada_attention_weights.R fit "$SLURM_ARRAY_TASK_ID" ${SMOKE:+--smoke}
echo "Exit $? at $(date)"
