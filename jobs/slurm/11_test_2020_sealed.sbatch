#!/bin/bash
#SBATCH --job-name=sar_phase3_2020
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=180G
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm-%j-sar_phase3_2020.out
#SBATCH --error=logs/slurm-%j-sar_phase3_2020.err

source /srv/software/easybuild/software/Anaconda3/2020.11/etc/profile.d/conda.sh
conda activate smap-xai-audit

cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

echo "=== Starting Phase 3 Sealed 2020 Test SLURM Job ==="
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-56}"
echo "Start time: $(date)"

python -u scripts/11_test_2020_sealed.py

echo "=== Phase 3 Sealed 2020 Test Finished Successfully at $(date) ==="
