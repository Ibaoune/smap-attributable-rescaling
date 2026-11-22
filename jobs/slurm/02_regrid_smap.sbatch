#!/bin/bash
#SBATCH --job-name=sar_regrid_smap
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=2016-2020%5
#SBATCH --output=logs/slurm-%A_%a-sar_regrid_smap.out
#SBATCH --error=logs/slurm-%A_%a-sar_regrid_smap.err

source /srv/software/easybuild/software/Anaconda3/2020.11/etc/profile.d/conda.sh
conda activate smap-xai-audit

cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

echo "Starting SMAP raw regridding for year ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date)"
python -u scripts/02_regrid_smap_raw.py --year ${SLURM_ARRAY_TASK_ID}
echo "Finished at $(date)"
