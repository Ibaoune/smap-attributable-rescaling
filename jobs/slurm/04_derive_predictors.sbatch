#!/bin/bash
#SBATCH --job-name=sar_predictors
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm-%j-sar_predictors.out
#SBATCH --error=logs/slurm-%j-sar_predictors.err

source /srv/software/easybuild/software/Anaconda3/2020.11/etc/profile.d/conda.sh
conda activate smap-xai-audit

cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

echo "Starting predictors derivation on $(hostname) at $(date)"
python -u scripts/04_derive_new_predictors.py
echo "Finished at $(date)"
