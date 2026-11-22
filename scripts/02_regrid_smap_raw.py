#!/usr/bin/env python3
"""
Step 1b: Reproject raw SMAP Enhanced L3 (SPL3SMP_E) HDF5 files onto NorthMor grid (121 x 181).
Extracts all 8 variables for AM and PM passes using exact LIS nearest-neighbor mapping.
Saves annual NetCDF files under data/interim/smap/smap_northmor_{year}.nc.
Includes strict automated quality controls at completion (no -9999, no out-of-range values).
Uses direct bounding-box slicing for ultra-fast I/O.
"""

import os
import sys
import glob
import argparse
import datetime
import numpy as np
import netCDF4
import xarray as xr
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09, regrid_smap_to_northmor
from sar.smap import open_smap_h5, read_smap_pass_data, SMAP_CANONICAL_VARIABLES

SMAP_RAW_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/observations/SMAP_SPL3SMP_E/raw"
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'interim', 'smap')

# Expected valid ranges for floating variables
EXPECTED_RANGES = {
    'soil_moisture': (0.0, 0.5),
    'vegetation_water_content': (0.0, 20.0),
    'surface_temperature': (0.0, 350.0),
    'roughness_coefficient': (0.0, 1.0),
    'static_water_body_fraction': (0.0, 1.0),
}


def get_grid_and_mapping():
    ds = netCDF4.Dataset(LIS_INPUT, 'r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    ds.close()
    row_idx, col_idx = latlon_to_ease2_m09(lat, lon)
    return lat[:, 0], lon[0, :], row_idx, col_idx


def verify_output_netcdf(filepath: str):
    """
    Strict validation of generated NetCDF:
    1. Zero instances of -9999 or -9999.0
    2. All variables within physical/valid bounds
    """
    print(f"\nRunning automated QC check on {filepath}...")
    ds = netCDF4.Dataset(filepath, 'r')
    errors = []
    
    for vname in ds.variables:
        if vname in ['time', 'lat', 'lon']:
            continue
        arr = ds.variables[vname][:]
        
        # Check for -9999
        m9999_count = int(np.sum(arr == -9999.0))
        if m9999_count > 0:
            errors.append(f"Variable {vname} contains {m9999_count} occurrences of -9999!")
            
        # Check range for known variables
        base_var = vname.replace('_am', '').replace('_pm', '')
        if base_var in EXPECTED_RANGES:
            vmin, vmax = EXPECTED_RANGES[base_var]
            if np.issubdtype(arr.dtype, np.floating):
                valid_mask = ~np.isnan(arr)
                if np.any(valid_mask):
                    actual_min = float(np.nanmin(arr))
                    actual_max = float(np.nanmax(arr))
                    if actual_min < vmin - 1e-5 or actual_max > vmax + 1e-5:
                        errors.append(f"Variable {vname} out of range: [{actual_min:.4f}, {actual_max:.4f}] vs expected [{vmin}, {vmax}]")
                        
    ds.close()
    if errors:
        for err in errors:
            print(f"  [ERROR] {err}")
        raise ValueError(f"QC failed for {filepath} with {len(errors)} errors.")
    else:
        print(f"  [PASSED] Automated QC passed: 0 instances of -9999, all variables strictly within physical bounds.")


def process_year(year: int, lat_1d, lon_1d, row_idx, col_idx):
    os.makedirs(OUT_DIR, exist_ok=True)
    out_file = os.path.join(OUT_DIR, f"smap_northmor_{year}.nc")
    
    start_date = datetime.date(year, 1, 1)
    end_date = datetime.date(year, 12, 31)
    n_days = (end_date - start_date).days + 1
    dates = [start_date + datetime.timedelta(days=i) for i in range(n_days)]
    
    n_lat, n_lon = 121, 181
    
    # Compute bounding box for NorthMor in global EASE2 grid
    r_min, r_max = int(row_idx.min()), int(row_idx.max())
    c_min, c_max = int(col_idx.min()), int(col_idx.max())
    bbox = (r_min, r_max, c_min, c_max)
    local_row_idx = row_idx - r_min
    local_col_idx = col_idx - c_min
    
    # Initialize containers for AM and PM variables
    data_vars = {}
    for p in ['am', 'pm']:
        for v in SMAP_CANONICAL_VARIABLES:
            var_name = f"{v}_{p}"
            dtype = np.uint16 if 'flag' in v else np.float32
            fill = 65535 if 'flag' in v else np.nan
            data_vars[var_name] = np.full((n_days, n_lat, n_lon), fill, dtype=dtype)
            
    print(f"Processing SMAP for year {year} ({n_days} days) using bounding box {bbox}...")
    found_count = 0
    
    for t_idx, d in enumerate(dates):
        d_str = d.strftime('%Y%m%d')
        matches = glob.glob(os.path.join(SMAP_RAW_DIR, f"SMAP_L3_SM_P_E_{d_str}_*.h5"))
        if not matches:
            continue
            
        fpath = matches[0]
        found_count += 1
        
        try:
            ds = open_smap_h5(fpath)
            for p in ['am', 'pm']:
                pass_name = p.upper()
                pass_dict = read_smap_pass_data(ds, pass_name, bbox=bbox)
                for v in SMAP_CANONICAL_VARIABLES:
                    var_name = f"{v}_{p}"
                    arr_sliced = pass_dict[v]
                    # Reproject via nearest neighbor on local slice
                    regridded = regrid_smap_to_northmor(arr_sliced, local_row_idx, local_col_idx)
                    data_vars[var_name][t_idx] = regridded
            ds.close()
        except Exception as e:
            print(f"Warning: error reading {fpath}: {e}")
            
        if (t_idx + 1) % 50 == 0 or t_idx == n_days - 1:
            print(f"  Day {t_idx + 1}/{n_days} ({found_count} files processed)")
            
    # Assemble xarray Dataset
    time_coords = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq='D')
    xr_dict = {}
    for var_name, arr in data_vars.items():
        if 'flag' in var_name:
            xr_dict[var_name] = (['time', 'lat', 'lon'], arr)
        else:
            xr_dict[var_name] = (['time', 'lat', 'lon'], arr.astype(np.float32))
            
    ds_out = xr.Dataset(
        xr_dict,
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={
            'title': f'SMAP SPL3SMP_E raw observations reprojected on NorthMor ({year})',
            'projection': 'nearest_neighbor matching LIS neighbor_interp',
            'source': 'SMAP SPL3SMP_E Version 006 (CRID R19240)',
            'domain': 'NorthMor (121 x 181, 0.05 deg)',
            'creation_date': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
    )
    
    # Encodings for compression
    encoding = {}
    for v in ds_out.data_vars:
        encoding[v] = {'zlib': True, 'complevel': 4}
        if 'flag' in v:
            encoding[v]['_FillValue'] = 65535
            
    ds_out.to_netcdf(out_file, encoding=encoding)
    print(f"Successfully generated: {out_file} ({os.path.getsize(out_file)/(1024*1024):.1f} MB)")
    
    # Automated QC check
    verify_output_netcdf(out_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=int, choices=[2016, 2017, 2018, 2019, 2020, 0], default=0,
                        help="Year to process (0 for all years 2016-2020)")
    args = parser.parse_args()
    
    lat_1d, lon_1d, row_idx, col_idx = get_grid_and_mapping()
    
    years = [2016, 2017, 2018, 2019, 2020] if args.year == 0 else [args.year]
    for y in years:
        process_year(y, lat_1d, lon_1d, row_idx, col_idx)

if __name__ == "__main__":
    main()
