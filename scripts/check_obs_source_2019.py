#!/usr/bin/env python3
"""
Verification of observation sources actually assimilated in 2019 DA experiments.
Compares DAOBS unscaled observations (Record 1) against:
1. Raw SMAP reprojected to 5 km.
2. V2 SMAP reprojected to 5 km.
Outputs full table to reports/check_obs_source_2019.csv.
No Decision 16 diagnostics computed.
"""

import os
import sys
import glob
import re
import netCDF4 as nc
import numpy as np

np.typeDict = np.sctypeDict
import h5py

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09, regrid_smap_to_northmor

WORKSPACE = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling"
P1_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco"
LIS_INPUT = os.path.join(P1_DIR, "data/lis_input/lis_input_NorthMor_5km.nc")
SMAP_RAW_DIR = os.path.join(P1_DIR, "data/observations/SMAP_SPL3SMP_E")
V2_OBS_DIR = os.path.join(WORKSPACE, "experiments/da_v2/obs")
REPORT_CSV = os.path.join(WORKSPACE, "reports/check_obs_source_2019.csv")

EXPERIMENTS = {
    "DA_nocdf_2019": "RAW",
    "DA_cdftrain_2019": "RAW",
    "DA_v2_2019": "V2",
}

def get_grid_mapping():
    ds = nc.Dataset(LIS_INPUT, "r")
    lat = ds["lat"][:]
    lon = ds["lon"][:]
    ds.close()
    return latlon_to_ease2_m09(lat, lon)

def read_daobs_rec1(filepath, n_lat=121, n_lon=181):
    try:
        with open(filepath, "rb") as f:
            hdr1 = np.fromfile(f, dtype=">u4", count=1)
            if len(hdr1) == 0:
                return None
            data = np.fromfile(f, dtype=">f4", count=n_lat * n_lon)
            if len(data) != n_lat * n_lon:
                return None
            return data.reshape((n_lat, n_lon))
    except Exception:
        return None

def main():
    print("Loading grid mapping...")
    row_map, col_map = get_grid_mapping()
    
    # Cache 5 km grids for raw and v2 by (date_str, pass_name)
    h5_cache_raw = {}
    h5_cache_v2 = {}
    
    def get_5km_grid(date_str, pass_name, source_type):
        cache = h5_cache_raw if source_type == "RAW" else h5_cache_v2
        key = (date_str, pass_name)
        if key in cache:
            return cache[key]
        
        yyyy = date_str[:4]
        mm = date_str[4:6]
        dd = date_str[6:8]
        date_dir = f"{yyyy}.{mm}.{dd}"
        
        if source_type == "RAW":
            patt = os.path.join(SMAP_RAW_DIR, date_dir, f"SMAP_L3_SM_P_E_{date_str}_*.h5")
        else:
            patt = os.path.join(V2_OBS_DIR, date_dir, f"SMAP_L3_SM_P_E_{date_str}_*.h5")
            
        matches = glob.glob(patt)
        if not matches:
            cache[key] = None
            return None
            
        fpath = matches[0]
        try:
            with h5py.File(fpath, "r") as hf:
                grp = f"Soil_Moisture_Retrieval_Data_{pass_name}"
                if grp not in hf:
                    cache[key] = None
                    return None
                ds_name = "soil_moisture" if pass_name == "AM" else "soil_moisture_pm"
                if ds_name not in hf[grp]:
                    cache[key] = None
                    return None
                sm_2d = hf[grp][ds_name][:]
                grid_5km = regrid_smap_to_northmor(sm_2d, row_map, col_map)
                cache[key] = grid_5km
                return grid_5km
        except Exception as e:
            print(f"Error reading {fpath}: {e}")
            cache[key] = None
            return None

    records = []
    
    for exp, exp_expected in EXPERIMENTS.items():
        print(f"Processing {exp}...")
        patt = os.path.join(WORKSPACE, "experiments/da_v2", exp, "output/**/DAOBS/**/LISDAOBS_*.1gs4r")
        files = sorted(glob.glob(patt, recursive=True))
        print(f"  Found {len(files)} DAOBS files.")
        
        for fpath in files:
            fname = os.path.basename(fpath)
            m = re.search(r"LISDAOBS_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})", fname)
            if not m:
                continue
            yyyy, mm, dd, hh, mn = m.groups()
            date_str = f"{yyyy}{mm}{dd}"
            hour = int(hh)
            pass_name = "AM" if hour < 12 else "PM"
            cycle_str = f"{yyyy}-{mm}-{dd} {hh}:{mn}"
            
            daobs_rec1 = read_daobs_rec1(fpath)
            if daobs_rec1 is None:
                continue
                
            raw_5km = get_5km_grid(date_str, pass_name, "RAW")
            v2_5km = get_5km_grid(date_str, pass_name, "V2")
            
            # Valid pixels where daobs has valid observation (0 < obs < 1)
            val_mask = (daobs_rec1 > 0) & (daobs_rec1 < 1.0)
            n_val = int(val_mask.sum())
            
            if n_val == 0:
                records.append({
                    "experiment": exp,
                    "cycle": cycle_str,
                    "pass": pass_name,
                    "expected_source": exp_expected,
                    "n_valid_pixels": 0,
                    "max_diff_raw": 0.0,
                    "max_diff_v2": 0.0,
                    "identified_source": "NO_OBS",
                    "status": "NO_OBS"
                })
                continue
                
            # Compare vs RAW
            if raw_5km is not None:
                diff_raw = np.abs(daobs_rec1[val_mask] - raw_5km[val_mask])
                max_diff_raw = float(diff_raw.max())
            else:
                max_diff_raw = 9999.0
                
            # Compare vs V2
            if v2_5km is not None:
                diff_v2 = np.abs(daobs_rec1[val_mask] - v2_5km[val_mask])
                max_diff_v2 = float(diff_v2.max())
            else:
                max_diff_v2 = 9999.0
                
            # Identify source
            is_raw = max_diff_raw < 1e-4
            is_v2 = max_diff_v2 < 1e-4
            
            if is_raw and is_v2:
                identified = "RAW_AND_V2_IDENTICAL"
                status = "CORRECT"
            elif is_raw and not is_v2:
                identified = "RAW"
                status = "CORRECT" if exp_expected == "RAW" else "WRONG_SOURCE"
            elif is_v2 and not is_raw:
                identified = "V2"
                status = "CORRECT" if exp_expected == "V2" else "WRONG_SOURCE"
            else:
                identified = "MISMATCH_OR_CORRUPT"
                status = "WRONG_SOURCE"
                
            records.append({
                "experiment": exp,
                "cycle": cycle_str,
                "pass": pass_name,
                "expected_source": exp_expected,
                "n_valid_pixels": n_val,
                "max_diff_raw": round(max_diff_raw, 6),
                "max_diff_v2": round(max_diff_v2, 6),
                "identified_source": identified,
                "status": status
            })

    # Write CSV
    os.makedirs(os.path.dirname(REPORT_CSV), exist_ok=True)
    header = "experiment,cycle,pass,expected_source,n_valid_pixels,max_diff_raw,max_diff_v2,identified_source,status\n"
    with open(REPORT_CSV, "w") as f:
        f.write(header)
        for r in records:
            line = f"{r['experiment']},{r['cycle']},{r['pass']},{r['expected_source']},{r['n_valid_pixels']},{r['max_diff_raw']},{r['max_diff_v2']},{r['identified_source']},{r['status']}\n"
            f.write(line)
            
    print(f"Report written to {REPORT_CSV}")
    
    # Print summary
    print("\n=== SUMMARY BY EXPERIMENT ===")
    import pandas as pd
    df = pd.DataFrame(records)
    for exp in EXPERIMENTS:
        df_exp = df[df['experiment'] == exp]
        total_cycles = len(df_exp)
        obs_cycles = len(df_exp[df_exp['n_valid_pixels'] > 0])
        correct = len(df_exp[df_exp['status'] == 'CORRECT'])
        wrong = len(df_exp[df_exp['status'] == 'WRONG_SOURCE'])
        print(f"Experiment {exp}:")
        print(f"  Total cycles: {total_cycles}, with valid obs: {obs_cycles}")
        print(f"  CORRECT: {correct}, WRONG_SOURCE: {wrong}")
        if wrong > 0:
            print("  Sample wrong source cycles:")
            print(df_exp[df_exp['status'] == 'WRONG_SOURCE'][['cycle', 'expected_source', 'identified_source', 'max_diff_raw', 'max_diff_v2']].head(10))

if __name__ == "__main__":
    main()
