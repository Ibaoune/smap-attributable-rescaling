#!/bin/bash
# Author: M. El Aabaribaoune (@um6p)
#SBATCH --job-name=DA_v2_2019_2019-01
#SBATCH --output=/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/logs/012019/DA_v2_2019_2019-01_%j.log
#SBATCH --error=/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/logs/012019/DA_v2_2019_2019-01_%j.err
#SBATCH --time=12:00:00
#SBATCH --account=empowermed-ahl6xm8o7mg-DEFAULT-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=32

set -e
cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco

source arch/arch_toubkal.env

echo "============================================="
echo " LIS NorthMor - DA - Period 2019-01"
echo " Start: 2019-01-01 | End: 2019-02-01"
echo "============================================="

time mpirun -n 32 ./src/lisf/lis/LIS -f /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/generated/configs/lis_DA_v2_2019_2019-01.config

echo "============================================="
echo " Post-processing & Validation "
echo "============================================="

SURF_RST=$(find /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/output/2019-01/SURFACEMODEL -name "LIS_RST_NOAHMP401_201902010000.d01.nc" | head -n 1)
PERT_RST=$(find /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/output/2019-01/DAPERT -name "LIS_DAPERT_201902010000.d01.bin" | head -n 1)

if [ -n "$SURF_RST" ] && [ -f "$SURF_RST" ] && [ -n "$PERT_RST" ] && [ -f "$PERT_RST" ]; then
    echo "Found restart files."
    cp "$SURF_RST" /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/restarts/surf/
    cp "$PERT_RST" /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/restarts/pert/
    echo "Successfully copied restarts to /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019/restarts/"
else
    echo "ERROR: Restart files not found in output directory!"
    exit 1
fi

echo "============================================="
echo " Submitting next period "
echo "============================================="
cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_v2_2019
echo "Chain complete!"
