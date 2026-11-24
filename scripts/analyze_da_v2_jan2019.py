#!/usr/bin/env python3
# Author: M. El Aabaribaoune (@um6p)
"""
Technical Control for DA_v2_2019 January 2019 Simulation.
Calculates strictly technical metrics day-by-day:
1. Expected 5 km pixels whose EASE cell contains a valid V2 observation
2. Non-missing innovations in EnKF files (innov_01 > -9900)
3. Ratio of assimilated innovations to expected pixels
4. Mean and standard deviation of innovations (O - F)
5. Layer 1 ensemble spread (mean over domain land pixels)
Strict rule: No interpretation of the V2 effect.
"""

import os
import sys
import glob
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import h5py
import netCDF4 as nc

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
OBS_DIR = os.path.join(BASE_DIR, "experiments", "da_v2", "obs")
EXP_DIR = os.path.join(BASE_DIR, "experiments", "da_v2", "DA_v2_2019")
ENKF_DIR = os.path.join(EXP_DIR, "output", "2019-01", "EnKF", "201901")
OUT_CSV = os.path.join(BASE_DIR, "reports", "check_da_v2_jan2019_technical.csv")

def main():
    print("=" * 80)
    print(" CONTRÔLE TECHNIQUE DA_v2_2019 - JANVIER 2019")
    print("=" * 80)

    # 1. Load domain coordinates and landmask
    with nc.Dataset(LIS_INPUT, 'r') as ds:
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
        landmask = ds.variables['LANDMASK'][:]
    
    row_map, col_map = latlon_to_ease2_m09(lat, lon)
    is_land = (landmask > 0.5)
    total_land_pixels = np.sum(is_land)
    print(f"Domaine NorthMor : {lat.shape[0]}x{lat.shape[1]} grilles, {total_land_pixels} pixels de terre valides.")

    records = []
    
    start_date = datetime(2019, 1, 1)
    end_date = datetime(2019, 1, 31)
    curr = start_date

    while curr <= end_date:
        date_iso = curr.strftime("%Y-%m-%d")
        yyyy = curr.strftime("%Y")
        mm = curr.strftime("%m")
        dd = curr.strftime("%d")

        # --- A. Expected 5 km observations from V2 HDF5 ---
        h5_path = os.path.join(OBS_DIR, f"{yyyy}.{mm}.{dd}", f"SMAP_L3_SM_P_E_{yyyy}{mm}{dd}_R19240_001.h5")
        
        n_exp_am = 0
        n_exp_pm = 0
        n_exp_total = 0
        n_exp_unique = 0
        
        if os.path.exists(h5_path):
            with h5py.File(h5_path, 'r') as hf:
                # AM Pass
                if 'Soil_Moisture_Retrieval_Data_AM' in hf:
                    sm_am = hf['Soil_Moisture_Retrieval_Data_AM/soil_moisture'][:]
                    q_am = hf['Soil_Moisture_Retrieval_Data_AM/retrieval_qual_flag'][:]
                    vwc_am = hf['Soil_Moisture_Retrieval_Data_AM/vegetation_water_content'][:]
                    valid_ease_am = (sm_am > -9900) & ((q_am & 1) == 0) & (vwc_am <= 5.0)
                    valid_5km_am = valid_ease_am[row_map, col_map] & is_land
                    n_exp_am = int(np.sum(valid_5km_am))
                else:
                    valid_5km_am = np.zeros_like(is_land, dtype=bool)

                # PM Pass
                if 'Soil_Moisture_Retrieval_Data_PM' in hf:
                    sm_pm = hf['Soil_Moisture_Retrieval_Data_PM/soil_moisture_pm'][:]
                    q_pm = hf['Soil_Moisture_Retrieval_Data_PM/retrieval_qual_flag_pm'][:]
                    vwc_pm = hf['Soil_Moisture_Retrieval_Data_PM/vegetation_water_content_pm'][:]
                    valid_ease_pm = (sm_pm > -9900) & ((q_pm & 1) == 0) & (vwc_pm <= 5.0)
                    valid_5km_pm = valid_ease_pm[row_map, col_map] & is_land
                    n_exp_pm = int(np.sum(valid_5km_pm))
                else:
                    valid_5km_pm = np.zeros_like(is_land, dtype=bool)

                n_exp_total = n_exp_am + n_exp_pm
                n_exp_unique = int(np.sum(valid_5km_am | valid_5km_pm))

        # --- B. Non-missing innovations from EnKF ---
        innov_pattern = os.path.join(ENKF_DIR, f"LIS_DA_EnKF_{yyyy}{mm}{dd}*_innov.a01.d01.nc")
        innov_files = sorted(glob.glob(innov_pattern))
        
        daily_innovations = []
        for inf in innov_files:
            try:
                with nc.Dataset(inf, 'r') as ids:
                    inv = ids.variables['innov_01'][:]
                    valid_inv = inv[inv > -9900]
                    if len(valid_inv) > 0:
                        daily_innovations.extend(valid_inv.tolist())
            except Exception as e:
                print(f"Erreur lecture innov {inf}: {e}")

        n_innov = len(daily_innovations)
        if n_innov > 0:
            arr_inv = np.array(daily_innovations, dtype=np.float64)
            innov_mean = float(np.mean(arr_inv))
            innov_std = float(np.std(arr_inv))
        else:
            innov_mean = np.nan
            innov_std = np.nan

        # Ratios
        ratio_total = (n_innov / n_exp_total) if n_exp_total > 0 else np.nan
        ratio_unique = (n_innov / n_exp_unique) if n_exp_unique > 0 else np.nan

        # --- C. Ensemble spread of Layer 1 ---
        # Spread is archived daily at 00:00:
        # LIS_DA_EnKF_{next_date}0000_spread.a01.d01.nc represents the state spread at 24:00 of current day
        next_dt = curr + timedelta(days=1)
        spread_matches = glob.glob(os.path.join(EXP_DIR, 'output', '2019-01', 'EnKF', '*', f'LIS_DA_EnKF_{next_dt.strftime("%Y%m%d")}0000_spread.a01.d01.nc'))
        spread_path = spread_matches[0] if spread_matches else ''
        
        spread_mean = np.nan
        spread_std = np.nan
        if os.path.exists(spread_path):
            try:
                with nc.Dataset(spread_path, 'r') as sds:
                    vname = 'ensspread_Soil Moisture Layer 1_01'
                    if vname in sds.variables:
                        sp = sds.variables[vname][:]
                        valid_sp = sp[(sp > -9900) & is_land]
                        if len(valid_sp) > 0:
                            spread_mean = float(np.mean(valid_sp))
                            spread_std = float(np.std(valid_sp))
            except Exception as e:
                print(f"Erreur lecture spread {spread_path}: {e}")

        # Optional: Start-of-day spread (at 00:00 of curr day)
        start_spread_path = os.path.join(ENKF_DIR, f"LIS_DA_EnKF_{yyyy}{mm}{dd}0000_spread.a01.d01.nc")
        start_spread_mean = np.nan
        if os.path.exists(start_spread_path):
            try:
                with nc.Dataset(start_spread_path, 'r') as sds:
                    vname = 'ensspread_Soil Moisture Layer 1_01'
                    if vname in sds.variables:
                        sp = sds.variables[vname][:]
                        valid_sp = sp[(sp > -9900) & is_land]
                        if len(valid_sp) > 0:
                            start_spread_mean = float(np.mean(valid_sp))
            except Exception:
                pass

        records.append({
            'date': date_iso,
            'expected_5km_am': n_exp_am,
            'expected_5km_pm': n_exp_pm,
            'expected_5km_total': n_exp_total,
            'expected_5km_unique': n_exp_unique,
            'n_innov': n_innov,
            'ratio_innov_vs_total': ratio_total,
            'ratio_innov_vs_unique': ratio_unique,
            'innov_mean': innov_mean,
            'innov_std': innov_std,
            'spread_l1_mean': spread_mean,
            'spread_l1_std': spread_std,
            'start_spread_l1_mean': start_spread_mean
        })

        curr += timedelta(days=1)

    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nResultats sauvegardes dans {OUT_CSV}")

    # Display preview
    print("\nAperçu des métriques calculées :")
    cols_disp = ['date', 'expected_5km_total', 'expected_5km_unique', 'n_innov', 'ratio_innov_vs_total', 'ratio_innov_vs_unique', 'innov_mean', 'innov_std', 'spread_l1_mean']
    print(df[cols_disp].to_string(index=False))

if __name__ == "__main__":
    main()
