#!/bin/bash
# Author: M. El Aabaribaoune (@um6p)
#SBATCH --job-name=DA_cdftrain_2019_2019-11
#SBATCH --output=/srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/logs/112019/DA_cdftrain_2019_2019-11_%j.log
#SBATCH --error=/srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/logs/112019/DA_cdftrain_2019_2019-11_%j.err
#SBATCH --time=12:00:00
#SBATCH --account=empowermed-ahl6xm8o7mg-DEFAULT-CPU
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=32

set -e
RUN_DIR="/srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/run"
mkdir -p "$RUN_DIR"
for d in data src arch shared domain; do
    if [ -d "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/$d" ]; then
        ln -sfn "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/$d" "$RUN_DIR/$d"
    fi
done
cd "$RUN_DIR"

source arch/arch_toubkal.env

echo "============================================="
echo " LIS NorthMor - DA - Period 2019-11"
echo " Start: 2019-11-01 | End: 2019-12-01"
echo "============================================="

time mpirun -n 32 ./src/lisf/lis/LIS -f /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/generated/configs/lis_DA_cdftrain_2019_2019-11.config

echo "============================================="
echo " Post-processing & Validation "
echo "============================================="

SURF_RST=$(find /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/output/2019-11/SURFACEMODEL -name "LIS_RST_NOAHMP401_201912010000.d01.nc" | head -n 1)
PERT_RST=$(find /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/output/2019-11/DAPERT -name "LIS_DAPERT_201912010000.d01.bin" | head -n 1)

if [ -n "$SURF_RST" ] && [ -f "$SURF_RST" ] && [ -n "$PERT_RST" ] && [ -f "$PERT_RST" ]; then
    echo "Found restart files."
    cp "$SURF_RST" /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/restarts/surf/
    cp "$PERT_RST" /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/restarts/pert/
    echo "Successfully copied restarts to /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019/restarts/"
else
    echo "ERROR: Restart files not found in output directory!"
    exit 1
fi

echo "============================================="
echo " Submitting next period "
echo "============================================="
cd /srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/experiments/da_v2/DA_cdftrain_2019
python3 scripts/chain_da.py --current-date 2019-12-01 --submit
