#!/usr/bin/env python3
# Author: M. El Aabaribaoune (@um6p)
"""
Control script for ENS_OL_2018Q4 restart at 2019-01-01:
- Checks both possible restart names (201901010000 and 201812312345)
- Computes Layer 1 domain-mean inter-member standard deviation (ddof=1)
- Computes Layer 1 domain-mean range (max - min)
- Compares with P1 DA restart (DA_nocdf_noirr_2016_2020) at 2019-01-01
Writes results to reports/check_ens_ol_2018q4.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import netCDF4 as nc

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
P1_RST = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/DA_nocdf_noirr_2016_2020/restarts/surf/LIS_RST_NOAHMP401_201901010000.d01.nc"
OUT_CSV = os.path.join(BASE_DIR, "reports/check_ens_ol_2018q4.csv")

# Potential restart paths for ENS_OL_2018Q4
RST_CANDIDATES = [
    os.path.join(BASE_DIR, "experiments/da_v2/ENS_OL_2018Q4/output/2018Q4/SURFACEMODEL/201901/LIS_RST_NOAHMP401_201901010000.d01.nc"),
    os.path.join(BASE_DIR, "experiments/da_v2/ENS_OL_2018Q4/output/2018Q4/SURFACEMODEL/201812/LIS_RST_NOAHMP401_201812312345.d01.nc"),
    os.path.join(BASE_DIR, "experiments/da_v2/ENS_OL_2018Q4/output/2018Q4/SURFACEMODEL/201810/LIS_RST_NOAHMP401_201810312345.d01.nc"),
]

def evaluate_restart(rst_path, name):
    if not os.path.exists(rst_path):
        return {
            'experiment': name,
            'restart_file': os.path.basename(rst_path),
            'status': 'NOT_FOUND',
            'std_l1_mean': np.nan,
            'range_l1_mean': np.nan
        }
    with nc.Dataset(rst_path) as ds:
        sm = ds.variables['SMC'][:]
        # sm has shape (4, 320600) -> layer 0 is layer 1 (20 members x 16030 tiles)
        sm_l1 = sm[0, :].reshape(20, 16030)
        std_members = np.std(sm_l1, axis=0, ddof=1)
        range_members = np.ptp(sm_l1, axis=0)
        return {
            'experiment': name,
            'restart_file': os.path.basename(rst_path),
            'status': 'OK',
            'std_l1_mean': float(np.mean(std_members)),
            'range_l1_mean': float(np.mean(range_members))
        }

def main():
    rows = []
    # 1. P1 restart
    rows.append(evaluate_restart(P1_RST, "P1_DA_nocdf_noirr"))
    
    # 2. ENS_OL_2018Q4 end restart
    found_final = False
    for cand in RST_CANDIDATES[:2]:
        if os.path.exists(cand):
            rows.append(evaluate_restart(cand, "ENS_OL_2018Q4_end"))
            found_final = True
            break
            
    if not found_final:
        # Check if intermediate 201810 restart exists
        if os.path.exists(RST_CANDIDATES[2]):
            rows.append(evaluate_restart(RST_CANDIDATES[2], "ENS_OL_2018Q4_month1_intermediate"))
        else:
            rows.append(evaluate_restart(RST_CANDIDATES[0], "ENS_OL_2018Q4_end"))

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(df.to_string(index=False))

if __name__ == '__main__':
    main()
