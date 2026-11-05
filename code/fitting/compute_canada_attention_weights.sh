#!/bin/bash
# SLURM array job for Rorqual (rorqual.alliancecan.ca): the ten brms
# attention-weight fits on the Bayesian model variants' simulated choices, one
# variant x range per task (see compute_canada_attention_weights.R). Each task
# takes one chiplet of a Rorqual node (8 cores: four chains, two threads each).
# Everything lives in one folder on the cluster, ~/links/scratch/attention_fits
# (Rorqual reaches scratch through the ~/links symlinks).
#
# 1. Laptop, from the project root:
#
#      rsync -av code/fitting/compute_canada_attention_weights.R code/fitting/compute_canada_attention_weights.sh results/compute_canada/inputs/*.rds beum@rorqual.alliancecan.ca:links/scratch/attention_fits/
#
# 2. Rorqual. The setup installs the R packages and CmdStan on the login node
#    (the compute nodes have no internet; it takes a while the first time and
#    is instant after), then the array is submitted:
#
#      cd ~/links/scratch/attention_fits && module load StdEnv/2023 r/4.4.0 && Rscript compute_canada_attention_weights.R setup && sbatch compute_canada_attention_weights.sh
#
#    squeue -u $USER shows the tasks; each writes log_<job>_<task>.txt here,
#    ending in "Done in ... h" and "Exit 0", and emails when it ends.
#
# 3. Laptop, when all ten are done (the 1.4 GB fits stay on the cluster):
#
#      rsync -av 'beum@rorqual.alliancecan.ca:links/scratch/attention_fits/outputs/*_draws.rds' results/compute_canada/outputs/
#      Rscript code/fitting/compute_canada_attention_weights.R collect
#
# Time: the same fit on the data took 6.6 hours on four M2 cores with one
# thread per chain. The 2.5-day limit is a wide safety margin (a task that
# hits it loses its fit); jobs up to 3 days get the same scheduling.

#SBATCH --account=def-webbr
#SBATCH --job-name=attn_sim
#SBATCH --array=1-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-12:00:00 # DD-HH:MM:SS; expected 5-10 h per task, this is a wide margin
#SBATCH --output=/home/beum/links/scratch/attention_fits/log_%A_%a.txt
#SBATCH --mail-user=b.eum@rotman.utoronto.ca
#SBATCH --mail-type=END,FAIL

# Go to the folder with the scripts and the inputs
cd /home/beum/links/scratch/attention_fits || exit 1

# Load modules
module load StdEnv/2023
module load r/4.4.0 2>/dev/null || module load r
export TMPDIR="$SLURM_TMPDIR"          # CmdStan's per-chain CSV output goes to the node's local SSD

echo "Task $SLURM_ARRAY_TASK_ID on $(hostname), $(date): $(R --version | head -1)"
Rscript compute_canada_attention_weights.R setup || exit 1     # only checks that the login-node setup is visible here
Rscript compute_canada_attention_weights.R fit "$SLURM_ARRAY_TASK_ID" ${SMOKE:+--smoke}
echo "Exit $? at $(date)"
