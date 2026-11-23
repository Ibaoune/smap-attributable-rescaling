#!/bin/bash
#SBATCH --job-name=sar_chk_saved_mod
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=180G
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=logs/slurm-%j-sar_chk_saved_mod.out
#SBATCH --error=logs/slurm-%j-sar_chk_saved_mod.err

source /srv/software/easybuild/software/Anaconda3/2020.11/etc/profile.d/conda.sh
conda activate smap-xai-audit

cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

echo "=== Starting Saved Model Check 2019 SLURM Job ==="
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-56}"
echo "Start time: $(date)"

python -u scripts/check_saved_model_2019.py

echo "=== Saved Model Check 2019 Finished Successfully at $(date) ==="
