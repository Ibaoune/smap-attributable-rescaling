#!/bin/bash
#SBATCH --job-name=sar_regrid_smap
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-%j-sar_regrid_smap.out
#SBATCH --error=logs/slurm-%j-sar_regrid_smap.err

source /srv/software/easybuild/software/Anaconda3/2020.11/etc/profile.d/conda.sh
conda activate smap-xai-audit

cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

echo "Starting SMAP raw regridding on $(hostname) at $(date)"
python scripts/02_regrid_smap_raw.py --year 0
echo "Finished at $(date)"
