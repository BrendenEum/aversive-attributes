#!/bin/bash
# SLURM array job for Fir (fir.alliancecan.ca): the ten brms attention-weight
# fits on the Bayesian model variants' simulated choices, one variant x range
# per task (see compute_canada_attention_weights.R). Each task takes one
# chiplet of a Fir node (8 cores: four chains, two threads each) and installs
# the R packages on its first run (Fir's compute nodes have internet; a lock
# makes one task install while the others wait), so there is nothing to set
# up by hand.
#
# 1. Laptop, from the project root (the inputs are not in git; -R keeps the
#    paths the R script relies on):
#
#      rsync -avR code/fitting/compute_canada_attention_weights.R code/fitting/compute_canada_attention_weights.sh results/compute_canada/inputs/ beum@fir.alliancecan.ca:scratch/aversive-attributes/
#
# 2. Fir:
#
#      cd ~/scratch/aversive-attributes && sbatch code/fitting/compute_canada_attention_weights.sh
#
#    If sbatch complains about the account, add --account=def-NAME, where
#    def-NAME is the folder in ~/project. squeue -u $USER shows the tasks;
#    each writes results/compute_canada/log_<job>_<task>.txt, ending in
#    "Done in ... h" and "Exit 0".
#
# 3. Laptop, when all ten are done (the 1.4 GB fits stay on Fir):
#
#      rsync -av 'beum@fir.alliancecan.ca:scratch/aversive-attributes/results/compute_canada/outputs/*_draws.rds' results/compute_canada/outputs/
#      Rscript code/fitting/compute_canada_attention_weights.R collect
#
# Time: the same fit on the data took 6.6 hours on four M2 cores with one
# thread per chain; the first task also spends up to an hour installing. The
# 24-hour limit is a safety margin (a task that hits it loses its fit).

#SBATCH --job-name=attn_sim
#SBATCH --array=1-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00
#SBATCH --output=results/compute_canada/log_%A_%a.txt

module load StdEnv/2023
module load r/4.4.0 2>/dev/null || module load r
export TMPDIR="$SLURM_TMPDIR"          # CmdStan's per-chain CSV output goes to the node's NVMe

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "Task $SLURM_ARRAY_TASK_ID on $(hostname), $(date): $(R --version | head -1)"
flock "$HOME/.compute_canada_setup.lock" Rscript code/fitting/compute_canada_attention_weights.R setup || exit 1
Rscript code/fitting/compute_canada_attention_weights.R fit "$SLURM_ARRAY_TASK_ID" ${SMOKE:+--smoke}
echo "Exit $? at $(date)"
