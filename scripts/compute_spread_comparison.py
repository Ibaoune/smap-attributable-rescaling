#!/usr/bin/env python3
# Author: M. El Aabaribaoune (@um6p)
"""
Spread Comparison Script (P1, DA_v2_2019, ENS_OPL_2015).
Computes exactly the metric spread_l1_mean on the 16,030 land pixels (landmask > 0.5):
- P1 DA_nocdf_noirr_2016_2020 for January 2019
- DA_v2_2019 for January 2019
- ENS_OPL_2015 last available spread reference
Outputs to reports/check_spread_comparison.csv
"""

import os
import glob
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import netCDF4 as nc

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
P1_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/DA_nocdf_noirr_2016_2020"
V2_DIR = os.path.join(BASE_DIR, "experiments", "da_v2", "DA_v2_2019")
OUT_CSV = os.path.join(BASE_DIR, "reports", "check_spread_comparison.csv")

def main():
    with nc.Dataset(LIS_INPUT) as ds:
        landmask = ds.variables['LANDMASK'][:]
    is_land = (landmask > 0.5)

    vname = 'ensspread_Soil Moisture Layer 1_01'

    # Last spread available in SPINUP (SPINUP_DA EnKF spread)
    spinup_last_file = '/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/SPINUP/SPINUP_DA/output/2015-12/EnKF/201601/LIS_DA_EnKF_201601010000_spread.a01.d01.nc'
    with nc.Dataset(spinup_last_file) as sds:
        sp = sds.variables[vname][:]
        spinup_spread_val = float(np.mean(sp[(sp > -9900) & is_land]))

    rows = []
    curr = datetime(2019, 1, 1)
    end = datetime(2019, 1, 31)

    while curr <= end:
        next_dt = curr + timedelta(days=1)
        dstr = curr.strftime('%Y-%m-%d')
        next_str = next_dt.strftime('%Y%m%d')
        
        # P1 spread
        p1_matches = glob.glob(os.path.join(P1_DIR, 'output', '2019-01', 'EnKF', '*', f'LIS_DA_EnKF_{next_str}0000_spread.a01.d01.nc'))
        p1_val = np.nan
        if p1_matches:
            with nc.Dataset(p1_matches[0]) as sds:
                sp = sds.variables[vname][:]
                p1_val = float(np.mean(sp[(sp > -9900) & is_land]))
                
        # DA_v2 spread
        v2_matches = glob.glob(os.path.join(V2_DIR, 'output', '2019-01', 'EnKF', '*', f'LIS_DA_EnKF_{next_str}0000_spread.a01.d01.nc'))
        v2_val = np.nan
        if v2_matches:
            with nc.Dataset(v2_matches[0]) as sds:
                sp = sds.variables[vname][:]
                v2_val = float(np.mean(sp[(sp > -9900) & is_land]))
                
        rows.append({
            'date': dstr,
            'range_l1_mean_P1': p1_val,
            'range_l1_mean_DA_v2_2019': v2_val,
            'range_l1_mean_SPINUP_DA_last': spinup_spread_val
        })
        curr += timedelta(days=1)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Generated {OUT_CSV}")

if __name__ == '__main__':
    main()
