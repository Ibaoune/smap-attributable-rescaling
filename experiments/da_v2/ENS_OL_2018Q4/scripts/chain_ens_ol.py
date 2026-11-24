#!/usr/bin/env python3
# Author: M. El Aabaribaoune (@um6p)
"""
Script to prepare and run the Ensemble Open Loop Spin-up for 2018 Q4 (ENS_OL_2018Q4)
to generate the 20-member spun-up ensemble restart at 2019-01-01.
"""

import argparse
import os
import sys
import subprocess
import configparser
from datetime import datetime

def read_config(config_path: str) -> configparser.SectionProxy:
    config = configparser.ConfigParser()
    if not os.path.exists(config_path):
        print(f"Error: Config file {config_path} not found.")
        sys.exit(1)
    config.read(config_path)
    return config['Experiment']

def main():
    parser = argparse.ArgumentParser(description="Prepare Ensemble OL for 2018 Q4.")
    parser.add_argument("--submit", action="store_true", help="Submit the job to SLURM immediately")
    args = parser.parse_args()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    os.chdir(base_dir)

    config_path = "config/experiment.ini"
    cfg = read_config(config_path)

    start_date = datetime.strptime(cfg['StartDate'], "%Y-%m-%d")
    end_date = datetime.strptime(cfg['EndDate'], "%Y-%m-%d")
    exp_name = cfg['ExperimentName']
    template_path = cfg['TemplateFile']
    start_mode = cfg.get('InitMode', 'restart')
    
    label = f"{start_date.strftime('%Y')}Q4"

    os.makedirs("generated/configs", exist_ok=True)
    os.makedirs("generated/jobs", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    os.makedirs(f"output/{label}", exist_ok=True)
    os.makedirs("restarts", exist_ok=True)
    
    # Check expanded 20-member restart file at 2018-10-01
    expected_start_restart = f"LIS_RST_NOAHMP401_{start_date.strftime('%Y%m%d')}0000_ens20.d01.nc"
    restart_path = os.path.join(base_dir, "restarts", expected_start_restart)
    if not os.path.exists(restart_path):
        alt_path = os.path.join(base_dir, "..", "restarts", expected_start_restart)
        if os.path.exists(alt_path):
            restart_path = os.path.abspath(alt_path)
        else:
            print(f"Error: Required spinup restart file {restart_path} not found!")
            sys.exit(1)

    with open(template_path, "r") as f:
        template_text = f.read()

    config_out = f"generated/configs/lis_{exp_name}.config"
    
    config_text = template_text
    config_text = config_text.replace("__START_MODE__", start_mode)
    config_text = config_text.replace("__START_YEAR__", start_date.strftime("%Y"))
    config_text = config_text.replace("__START_MONTH__", start_date.strftime("%m"))
    config_text = config_text.replace("__START_DAY__", start_date.strftime("%d"))
    config_text = config_text.replace("__END_YEAR__", end_date.strftime("%Y"))
    config_text = config_text.replace("__END_MONTH__", end_date.strftime("%m"))
    config_text = config_text.replace("__END_DAY__", end_date.strftime("%d"))
    config_text = config_text.replace("__SURF_RESTART_FILE__", restart_path)
    config_text = config_text.replace("__OUTPUT_DIR__", f"{base_dir}/output/{label}")
    config_text = config_text.replace("__DIAGNOSTIC_OUT__", f"{base_dir}/logs/lislog_{label}")

    with open(config_out, "w") as f:
        f.write(config_text)

    # Also symlink lis_ens_ol_2018Q4.config if needed
    alt_config_out = f"generated/configs/lis_ens_ol_2018Q4.config"
    with open(alt_config_out, "w") as f:
        f.write(config_text)

    # Expected Outputs at 2019-01-01
    expected_end_surf_restart = f"LIS_RST_NOAHMP401_{end_date.strftime('%Y%m%d')}0000.d01.nc"
    expected_end_pert_restart = f"LIS_DAPERT_{end_date.strftime('%Y%m%d')}0000.d01.bin"

    job_out = f"generated/jobs/job_{exp_name}.sh"

    lis_root = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco"

    slurm_script = f"""#!/bin/bash
# Author: M. El Aabaribaoune (@um6p)
#SBATCH --job-name={exp_name}
#SBATCH --output={base_dir}/logs/{exp_name}_%j.log
#SBATCH --error={base_dir}/logs/{exp_name}_%j.err
#SBATCH --time=12:00:00
#SBATCH --account=empowermed-ahl6xm8o7mg-DEFAULT-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=32

set -e
cd {lis_root}

source arch/arch_toubkal.env

echo "============================================="
echo " LIS NorthMor - {exp_name} - Period {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}"
echo "============================================="

time mpirun -n 32 ./src/lisf/lis/LIS -f {base_dir}/generated/configs/lis_{exp_name}.config

echo "============================================="
echo " Spin-up complete: 2019-01-01 restart ready"
echo "============================================="
"""

    with open(job_out, "w") as f:
        f.write(slurm_script)
    os.chmod(job_out, 0o755)

    print(f"Generated {config_out} and {job_out}")

    if args.submit:
        print(f"Submitting {job_out} to SLURM...")
        subprocess.run(["sbatch", job_out], check=True)

if __name__ == "__main__":
    main()
