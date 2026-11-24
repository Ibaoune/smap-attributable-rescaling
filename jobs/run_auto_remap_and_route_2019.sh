#!/bin/bash
#SBATCH --job-name=auto_route_2019
#SBATCH --output=logs/auto_route_2019_%j.log
#SBATCH --error=logs/auto_route_2019_%j.err
#SBATCH --time=36:00:00
#SBATCH --account=empowermed-ahl6xm8o7mg-DEFAULT-CPU
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=16

set -e
cd /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling

source /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/arch/arch_toubkal.env

python3 scripts/auto_remap_and_route_2019.py
