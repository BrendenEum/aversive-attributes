#!/bin/bash
# SLURM array job: the ten brms attention-weight fits on the Bayesian model
# variants' simulated choices (see compute_canada_attention_weights.R), one
# variant x range per task, four cores each (one chain per core).
#
# One-time setup on the cluster (R packages go to your home library; building
# CmdStan takes ten to twenty minutes, so do it in a short interactive job
# rather than on the login node):
#
#   module load StdEnv/2023 r/4.4.0        # `module spider r` lists the versions
#   salloc --time=1:00:00 --cpus-per-task=4 --mem=8G --account=def-CHANGEME
#   Rscript -e 'install.packages(c("brms", "posterior", "dplyr"), repos = "https://cloud.r-project.org")'
#   Rscript -e 'install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))'
#   Rscript -e 'cmdstanr::install_cmdstan(cores = 4)'
#   exit
#
# Getting the project there (code through git, inputs by copy since they are
# not tracked), from your laptop:
#
#   ssh USER@cedar.alliancecan.ca 'cd ~/scratch && git clone https://github.com/BrendenEum/aversive-attributes.git'
#   rsync -av results/compute_canada/inputs/ USER@cedar.alliancecan.ca:~/scratch/aversive-attributes/results/compute_canada/inputs/
#
# Submitting, from the project root on the cluster (fill in the account):
#
#   mkdir -p results/logs
#   sbatch code/fitting/compute_canada_attention_weights.sh
#   squeue -u $USER                          # progress; logs in results/logs/
#
# Bringing the compact files back (the full fits stay on the cluster):
#
#   rsync -av 'USER@cedar.alliancecan.ca:~/scratch/aversive-attributes/results/compute_canada/outputs/*_draws.rds' results/compute_canada/outputs/
#   Rscript code/fitting/compute_canada_attention_weights.R collect
#
# Time: the same fit on the data took 6.6 hours on four M2 cores; cluster
# cores are slower per core, so the limit below allows for twice that. Jobs
# under 24 hours schedule faster than longer ones.

#SBATCH --job-name=attn_sim
#SBATCH --account=def-CHANGEME
#SBATCH --array=1-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-20:00
#SBATCH --output=results/logs/compute_canada_%A_%a.log

module load StdEnv/2023 r/4.4.0

export TMPDIR="$SLURM_TMPDIR"          # CmdStan's per-chain CSV output goes to node-local disk
export OMP_NUM_THREADS=1               # one chain per core, no nested threading

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "Task $SLURM_ARRAY_TASK_ID on $(hostname), $(date)"
Rscript code/fitting/compute_canada_attention_weights.R fit "$SLURM_ARRAY_TASK_ID"
echo "Exit $? at $(date)"
