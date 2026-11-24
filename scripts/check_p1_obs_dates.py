#!/usr/bin/env python3
"""
Exhaustive verification of observation dates in P1 runs (2016-2020):
- DA_nocdf_noirr_2016_2020
- DA_cdf_noirr_2016_2020
Compares DAOBS unscaled observation (Record 1) to Raw SMAP reprojected to 5 km of the same date.
If any mismatch is found, tests adjacent dates.
Outputs full report to reports/check_p1_obs_dates.csv.
"""

import os
import sys
import glob
import re
import datetime
import netCDF4 as nc
import numpy as np

np.typeDict = np.sctypeDict
import h5py

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09, regrid_smap_to_northmor

P1_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco"
LIS_INPUT = os.path.join(P1_DIR, "data/lis_input/lis_input_NorthMor_5km.nc")
SMAP_RAW_DIR = os.path.join(P1_DIR, "data/observations/SMAP_SPL3SMP_E")
REPORT_CSV = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), "reports/check_p1_obs_dates.csv")

RUNS = {
    "DA_nocdf_noirr_2016_2020": os.path.join(P1_DIR, "experiments/NorthMor/matrix_2016_2020/DA_nocdf_noirr_2016_2020/output"),
    "DA_cdf_noirr_2016_2020": os.path.join(P1_DIR, "experiments/NorthMor/matrix_2016_2020/DA_cdf_noirr_2016_2020/output"),
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
    
    # Cache raw 5 km grids by (date_str, pass_name)
    raw_cache = {}
    
    def get_raw_5km(date_str, pass_name):
        key = (date_str, pass_name)
        if key in raw_cache:
            return raw_cache[key]
        yyyy = date_str[:4]
        mm = date_str[4:6]
        dd = date_str[6:8]
        patt = os.path.join(SMAP_RAW_DIR, f"{yyyy}.{mm}.{dd}", f"SMAP_L3_SM_P_E_{date_str}_*.h5")
        matches = glob.glob(patt)
        if not matches:
            raw_cache[key] = None
            return None
        fpath = matches[0]
        try:
            with h5py.File(fpath, "r") as hf:
                grp = f"Soil_Moisture_Retrieval_Data_{pass_name}"
                if grp not in hf:
                    raw_cache[key] = None
                    return None
                ds_name = "soil_moisture" if pass_name == "AM" else "soil_moisture_pm"
                if ds_name not in hf[grp]:
                    raw_cache[key] = None
                    return None
                sm_2d = hf[grp][ds_name][:]
                grid_5km = regrid_smap_to_northmor(sm_2d, row_map, col_map)
                raw_cache[key] = grid_5km
                return grid_5km
        except Exception:
            raw_cache[key] = None
            return None

    records = []
    
    for run_name, out_root in RUNS.items():
        print(f"\nProcessing {run_name}...")
        patt = os.path.join(out_root, "**/DAOBS/**/LISDAOBS_*.1gs4r")
        files = sorted(glob.glob(patt, recursive=True))
        print(f"  Found {len(files)} DAOBS files.")
        
        mismatch_count = 0
        match_count = 0
        no_obs_count = 0
        
        for idx, fpath in enumerate(files):
            if (idx + 1) % 500 == 0:
                print(f"  Checked {idx + 1}/{len(files)} files...")
                
            fname = os.path.basename(fpath)
            m = re.search(r"LISDAOBS_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})", fname)
            if not m:
                continue
            yyyy, mm, dd, hh, mn = m.groups()
            date_str = f"{yyyy}{mm}{dd}"
            hour = int(hh)
            pass_name = "AM" if hour < 12 else "PM"
            cycle_str = f"{yyyy}-{mm}-{dd} {hh}:{mn}"
            
            daobs = read_daobs_rec1(fpath)
            if daobs is None:
                continue
                
            val_mask = (daobs > 0) & (daobs < 1.0)
            n_val = int(val_mask.sum())
            
            if n_val == 0:
                no_obs_count += 1
                records.append({
                    "experiment": run_name,
                    "cycle": cycle_str,
                    "date": date_str,
                    "pass": pass_name,
                    "n_valid_pixels": 0,
                    "max_diff_same_date": 0.0,
                    "status": "NO_OBS",
                    "matched_date": date_str
                })
                continue
                
            raw_same = get_raw_5km(date_str, pass_name)
            if raw_same is not None:
                diff = np.abs(daobs[val_mask] - raw_same[val_mask])
                max_diff = float(diff.max())
            else:
                max_diff = 9999.0
                
            if max_diff < 1e-4:
                match_count += 1
                records.append({
                    "experiment": run_name,
                    "cycle": cycle_str,
                    "date": date_str,
                    "pass": pass_name,
                    "n_valid_pixels": n_val,
                    "max_diff_same_date": round(max_diff, 6),
                    "status": "MATCH_SAME_DATE",
                    "matched_date": date_str
                })
            else:
                mismatch_count += 1
                # Try surrounding dates (-3 to +3 days)
                cur_dt = datetime.datetime.strptime(date_str, "%Y%m%d")
                alt_matched = None
                for offset in [-2, -1, 1, 2, -3, 3]:
                    alt_dt = cur_dt + datetime.timedelta(days=offset)
                    alt_date_str = alt_dt.strftime("%Y%m%d")
                    raw_alt = get_raw_5km(alt_date_str, pass_name)
                    if raw_alt is not None:
                        diff_alt = np.abs(daobs[val_mask] - raw_alt[val_mask])
                        if float(diff_alt.max()) < 1e-4:
                            alt_matched = alt_date_str
                            break
                            
                status_str = f"MISMATCH_MATCHES_{alt_matched}" if alt_matched else "MISMATCH_NO_DATE_MATCH"
                print(f"  [MISMATCH] {run_name} {cycle_str} (n={n_val}, max_diff={max_diff:.4f}) -> {status_str}")
                records.append({
                    "experiment": run_name,
                    "cycle": cycle_str,
                    "date": date_str,
                    "pass": pass_name,
                    "n_valid_pixels": n_val,
                    "max_diff_same_date": round(max_diff, 6),
                    "status": status_str,
                    "matched_date": alt_matched if alt_matched else "NONE"
                })

        print(f"Finished {run_name}: MATCH={match_count}, NO_OBS={no_obs_count}, MISMATCH={mismatch_count}")

    # Write report
    os.makedirs(os.path.dirname(REPORT_CSV), exist_ok=True)
    with open(REPORT_CSV, "w") as f:
        f.write("experiment,cycle,date,pass,n_valid_pixels,max_diff_same_date,status,matched_date\n")
        for r in records:
            f.write(f"{r['experiment']},{r['cycle']},{r['date']},{r['pass']},{r['n_valid_pixels']},{r['max_diff_same_date']},{r['status']},{r['matched_date']}\n")
            
    print(f"\nAll results saved to {REPORT_CSV}")
    
    # Global summary
    import pandas as pd
    df = pd.DataFrame(records)
    print("\n=== GLOBAL SYNTHESIS ===")
    for run_name in RUNS:
        df_r = df[df['experiment'] == run_name]
        total = len(df_r)
        with_obs = len(df_r[df_r['n_valid_pixels'] > 0])
        matches = len(df_r[df_r['status'] == 'MATCH_SAME_DATE'])
        mismatches = len(df_r[df_r['status'].str.startswith('MISMATCH')])
        print(f"{run_name}:")
        print(f"  Total cycles: {total}")
        print(f"  Cycles with valid observations: {with_obs}")
        print(f"  Exact matches with same-date raw SMAP: {matches} ({matches/with_obs*100.0:.2f}%)")
        print(f"  Date mismatches / shifted observations: {mismatches}")

if __name__ == "__main__":
    main()
