#!/bin/bash
# Author: M. El Aabaribaoune (@um6p)
#SBATCH --job-name=ENS_OL_2018Q4
#SBATCH --output=/srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/ENS_OL_2018Q4/logs/ENS_OL_2018Q4_%j.log
#SBATCH --error=/srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/ENS_OL_2018Q4/logs/ENS_OL_2018Q4_%j.err
#SBATCH --time=12:00:00
#SBATCH --account=empowermed-ahl6xm8o7mg-DEFAULT-CPU
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=32

set -e
cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco

source arch/arch_toubkal.env

echo "============================================="
echo " LIS NorthMor - ENS_OL_2018Q4 - Period 2018-10-01 to 2019-01-01"
echo "============================================="

time mpirun -n 32 ./src/lisf/lis/LIS -f /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/ENS_OL_2018Q4/generated/configs/lis_ENS_OL_2018Q4.config

echo "============================================="
echo " Spin-up complete: 2019-01-01 restart ready"
echo "============================================="
