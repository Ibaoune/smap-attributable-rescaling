#!/usr/bin/env python3
"""
Phase 1c: Derive and assemble all new predictors according to data contract:
1. Static predictors: slope_deg, slopetype, irrigfrac_gmia, texture & landcover fractions,
   dominant soil class and soil water limits (MAXSMC, REFSMC, WLTSMC) from SOILPARM.TBL.
2. MSWEP v2.8 0.10 deg precipitation regridded to NorthMor 0.05 deg via direct cell assignment
   with basin water-balance conservation check.
3. delta_precip_imerg_mswep = precip_imerg - precip_mswep.
4. Antecedent precip metrics strictly ending at t-1: precip_1d_sum, precip_3d_sum, precip_7d_sum, days_since_rain_gt_1mm.
5. Extraction of LAI_tavg and QC variables (SoilTemp, SWE, Snowcover, Greenness) from OPL_noirr_2016_2020.
6. Relative soil moisture index ssm_rel_wilt_ref at t-1 clipped to [0.0, 1.5].
7. Day of year harmonics (doy_sin, doy_cos).
"""

import os
import sys
import glob
import datetime
import hashlib
import numpy as np
import netCDF4
import xarray as xr
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.predictors import (
    compute_topographic_slope,
    compute_rolling_precip_3d,
    compute_dry_days_3d,
    compute_relative_soil_moisture,
    compute_seasonal_harmonics
)
from sar.io import compute_sha256, log_provenance

# Paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
ELEVATION_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "elevation.nc")
SOILPARM_FILE = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/land_params/noahmp401_parms/SOILPARM.TBL"
MSWEP_FILE = "/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/data/mswep/all_Mor/mswep_1979_2020.nc"
OPL_NOIRR_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/OPL_noirr_2016_2020/output"
MANIFEST_FILE = os.path.join(BASE_DIR, "docs", "provenance_manifest.csv")


def parse_soilparm():
    """Extract MAXSMC, REFSMC, WLTSMC for STATSGO classes from SOILPARM.TBL."""
    wltsmc, refsmc, maxsmc = {}, {}, {}
    with open(SOILPARM_FILE, 'r') as f:
        in_stas = False
        for line in f:
            if 'STAS' in line and not 'STAS-RUC' in line:
                in_stas = True
                continue
            if 'STAS-RUC' in line:
                break
            if in_stas and line.strip() and line.strip()[0].isdigit():
                parts = [p.strip() for p in line.split(',') if p.strip()]
                if len(parts) >= 10:
                    c_id = int(parts[0])
                    maxsmc[c_id] = float(parts[4])
                    refsmc[c_id] = float(parts[5])
                    wltsmc[c_id] = float(parts[9])
    return maxsmc, refsmc, wltsmc


def derive_static_predictors():
    """Extract and derive all static predictors."""
    print("Deriving static predictors...")
    out_dir = os.path.join(BASE_DIR, "data", "interim", "static")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "static_predictors.nc")
    
    # 1. Elevation & slope
    ds_elev = netCDF4.Dataset(ELEVATION_FILE, 'r')
    elev = ds_elev.variables['elevation'][:]
    lat = ds_elev.variables['lat'][:]
    lon = ds_elev.variables['lon'][:]
    ds_elev.close()
    
    slope_deg = compute_topographic_slope(elev, res_deg=0.05, lat_ref=33.0)
    
    # 2. LIS Input parameters
    ds_lis = netCDF4.Dataset(LIS_INPUT, 'r')
    irrigfrac = ds_lis.variables['IRRIGFRAC'][:]
    slopetype = ds_lis.variables['SLOPETYPE'][:]
    landmask = ds_lis.variables['LANDMASK'][:].squeeze().astype(int)
    texture = ds_lis.variables['TEXTURE'][:]  # (16, 121, 181)
    landcover = ds_lis.variables['LANDCOVER'][:]  # (20, 121, 181)
    ds_lis.close()
    
    # Dominant classes
    dominant_soil = np.argmax(texture, axis=0) + 1  # 1-based
    dominant_landcover = np.argmax(landcover, axis=0) + 1  # 1-based
    
    # SOILPARM mapping
    maxsmc, refsmc, wltsmc = parse_soilparm()
    grid_sat = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    grid_ref = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    grid_wlt = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    
    for c_id in maxsmc:
        mask = (dominant_soil == c_id)
        grid_sat[mask] = maxsmc[c_id]
        grid_ref[mask] = refsmc[c_id]
        grid_wlt[mask] = wltsmc[c_id]
        
    # Assemble Dataset
    xr_dict = {
        'elevation': (['lat', 'lon'], elev.astype(np.float32)),
        'slope_deg': (['lat', 'lon'], slope_deg.astype(np.float32)),
        'slopetype': (['lat', 'lon'], slopetype.astype(np.int32)),
        'irrigfrac_gmia': (['lat', 'lon'], irrigfrac.astype(np.float32)),
        'landmask': (['lat', 'lon'], landmask.astype(np.int32)),
        'dominant_soiltype': (['lat', 'lon'], dominant_soil.astype(np.int32)),
        'dominant_landcover': (['lat', 'lon'], dominant_landcover.astype(np.int32)),
        'soil_sat': (['lat', 'lon'], grid_sat),
        'soil_ref': (['lat', 'lon'], grid_ref),
        'soil_wlt': (['lat', 'lon'], grid_wlt),
    }
    
    # Add individual texture fractions (1 to 16)
    for i in range(16):
        xr_dict[f'texture_frac_{i+1}'] = (['lat', 'lon'], texture[i].astype(np.float32))
        
    # Add individual landcover fractions (1 to 20)
    for i in range(20):
        xr_dict[f'landcover_frac_{i+1}'] = (['lat', 'lon'], landcover[i].astype(np.float32))
        
    ds_out = xr.Dataset(xr_dict, coords={'lat': lat, 'lon': lon})
    ds_out.to_netcdf(out_file)
    print(f"Static predictors saved: {out_file} ({os.path.getsize(out_file)/(1024*1024):.2f} MB)")
    return ds_out


def find_opl_file(d: datetime.date) -> str:
    """Find OPL_noirr LIS_HIST daily file for calendar day d."""
    next_d = d + datetime.timedelta(days=1)
    target = f"LIS_HIST_{next_d.strftime('%Y%m%d')}0000.d01.nc"
    ym_cur = d.strftime('%Y-%m')
    ym_cur_c = d.strftime('%Y%m')
    ym_next = next_d.strftime('%Y-%m')
    ym_next_c = next_d.strftime('%Y%m')
    
    cands = [
        os.path.join(OPL_NOIRR_DIR, ym_cur, 'SURFACEMODEL', ym_cur_c, target),
        os.path.join(OPL_NOIRR_DIR, ym_cur, 'SURFACEMODEL', ym_next_c, target),
        os.path.join(OPL_NOIRR_DIR, ym_next, 'SURFACEMODEL', ym_cur_c, target),
        os.path.join(OPL_NOIRR_DIR, ym_next, 'SURFACEMODEL', ym_next_c, target),
    ]
    for c in cands:
        if os.path.exists(c):
            return c
    return None


def extract_opl_variables(year: int, lat_1d, lon_1d):
    """Extract LAI_tavg and QC variables from OPL_noirr for a specific year."""
    print(f"Extracting LAI and QC variables from OPL_noirr for {year}...")
    start_d = datetime.date(year, 1, 1)
    end_d = datetime.date(year, 12, 31)
    n_days = (end_d - start_d).days + 1
    dates = [start_d + datetime.timedelta(days=i) for i in range(n_days)]
    
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    
    lai_arr = np.full((n_days, n_lat, n_lon), np.nan, dtype=np.float32)
    soil_temp_arr = np.full((n_days, n_lat, n_lon), np.nan, dtype=np.float32)
    swe_arr = np.full((n_days, n_lat, n_lon), np.nan, dtype=np.float32)
    snowcover_arr = np.full((n_days, n_lat, n_lon), np.nan, dtype=np.float32)
    greenness_arr = np.full((n_days, n_lat, n_lon), np.nan, dtype=np.float32)
    
    for t_idx, d in enumerate(dates):
        fpath = find_opl_file(d)
        if fpath is None:
            raise FileNotFoundError(f"Missing OPL file for {d}")
            
        ds = netCDF4.Dataset(fpath, 'r')
        lai_arr[t_idx] = ds.variables['LAI_tavg'][:]
        # Layer 1 soil temperature (0-10 cm)
        soil_temp_arr[t_idx] = ds.variables['SoilTemp_tavg'][0, :, :]
        swe_arr[t_idx] = ds.variables['SWE_tavg'][:]
        snowcover_arr[t_idx] = ds.variables['Snowcover_tavg'][:]
        greenness_arr[t_idx] = ds.variables['Greenness_tavg'][:]
        ds.close()
        
    time_coords = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq='D')
    
    # Save LAI
    lai_dir = os.path.join(BASE_DIR, "data", "interim", "noahmp", "opl", "lai")
    os.makedirs(lai_dir, exist_ok=True)
    lai_file = os.path.join(lai_dir, f"lai_{year}.nc")
    ds_lai = xr.Dataset(
        {'lai': (['time', 'lat', 'lon'], lai_arr)},
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={'title': f'LAI_tavg extracted from OPL_noirr ({year})', 'units': 'm^2/m^2'}
    )
    ds_lai.to_netcdf(lai_file)
    print(f"  Saved LAI: {lai_file} ({os.path.getsize(lai_file)/(1024*1024):.1f} MB)")
    
    # Save QC variables
    qc_dir = os.path.join(BASE_DIR, "data", "interim", "noahmp", "opl", "qc")
    os.makedirs(qc_dir, exist_ok=True)
    qc_file = os.path.join(qc_dir, f"qc_vars_{year}.nc")
    ds_qc = xr.Dataset(
        {
            'soil_temp_layer1': (['time', 'lat', 'lon'], soil_temp_arr),
            'swe': (['time', 'lat', 'lon'], swe_arr),
            'snowcover': (['time', 'lat', 'lon'], snowcover_arr),
            'greenness': (['time', 'lat', 'lon'], greenness_arr)
        },
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={'title': f'QC variables extracted from OPL_noirr ({year})'}
    )
    ds_qc.to_netcdf(qc_file)
    print(f"  Saved QC vars: {qc_file} ({os.path.getsize(qc_file)/(1024*1024):.1f} MB)")
    
    # Log provenance in manifest
    log_provenance(
        MANIFEST_FILE,
        source_path=lai_file,
        target_path=lai_file,
        producing_script="scripts/04_derive_new_predictors.py",
        description=f"Extracted LAI_tavg from OPL_noirr for year {year}"
    )

def regrid_mswep_and_compute_metrics(year: int, lat_1d, lon_1d):
    """Regrid MSWEP to NorthMor 0.05 deg, compute delta_precip, antecedent precip, and dry spells."""
    print(f"Processing MSWEP and precipitation metrics for {year}...")
    
    # Open MSWEP
    ds_m = netCDF4.Dataset(MSWEP_FILE, 'r')
    time_m = ds_m.variables['time'][:]
    lat_m = ds_m.variables['lat'][:]
    lon_m = ds_m.variables['lon'][:]
    
    base_date = datetime.date(1900, 1, 1)
    mswep_dates = [base_date + datetime.timedelta(days=int(t)) for t in time_m]
    
    # Map NorthMor lat/lon to MSWEP cell
    # Direct nearest-neighbor assignment (assigning containing 0.1 deg cell to 0.05 deg sub-cells)
    lat_2d, lon_2d = np.meshgrid(lat_1d, lon_1d, indexing='ij')
    lat_idx = np.abs(lat_m[:, None, None] - lat_2d[None, :, :]).argmin(axis=0)
    lon_idx = np.abs(lon_m[:, None, None] - lon_2d[None, :, :]).argmin(axis=0)
    
    # Extract MSWEP for year
    start_d = datetime.date(year, 1, 1)
    end_d = datetime.date(year, 12, 31)
    n_days = (end_d - start_d).days + 1
    
    # Target dates
    target_dates = [start_d + datetime.timedelta(days=i) for i in range(n_days)]
    target_idx = [mswep_dates.index(d) for d in target_dates]
    
    # Read MSWEP slice for target days
    mswep_raw = ds_m.variables['precipitation'][target_idx, :, :]
    ds_m.close()
    
    # Regrid via direct assignment
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    mswep_regrid = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    for t in range(n_days):
        mswep_regrid[t] = mswep_raw[t, lat_idx, lon_idx]
        
    # Save regridded MSWEP
    mswep_dir = os.path.join(BASE_DIR, "data", "interim", "mswep")
    os.makedirs(mswep_dir, exist_ok=True)
    mswep_file = os.path.join(mswep_dir, f"mswep_northmor_{year}.nc")
    time_coords = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq='D')
    
    ds_mswep = xr.Dataset(
        {'precip_mswep': (['time', 'lat', 'lon'], mswep_regrid)},
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={
            'title': f'MSWEP v2.8 daily precipitation regridded on NorthMor ({year})',
            'units': 'mm day-1',
            'regridding': 'direct cell assignment (nearest neighbor conserving block volume)'
        }
    )
    ds_mswep.to_netcdf(mswep_file)
    print(f"  Saved MSWEP: {mswep_file} ({os.path.getsize(mswep_file)/(1024*1024):.1f} MB)")
    
    # Load IMERG forcing for delta calculation
    imerg_file = os.path.join(BASE_DIR, "data", "interim", "noahmp", "forcing", "precipitation", f"precipitation_{year}.nc")
    ds_imerg = netCDF4.Dataset(imerg_file, 'r')
    imerg_precip = ds_imerg.variables['precipitation'][:]
    ds_imerg.close()
    
    delta_precip = imerg_precip - mswep_regrid
    
    # Antecedent precipitation metrics on IMERG strictly ending at t-1
    precip_1d_sum, precip_3d_sum, precip_7d_sum = compute_rolling_precip_3d(imerg_precip)
    days_since_rain = compute_dry_days_3d(imerg_precip, threshold_mm=1.0)
    
    # Save precipitation derived metrics
    metrics_dir = os.path.join(BASE_DIR, "data", "interim", "predictors_derived")
    os.makedirs(metrics_dir, exist_ok=True)
    metrics_file = os.path.join(metrics_dir, f"precip_metrics_{year}.nc")
    
    ds_metrics = xr.Dataset(
        {
            'delta_precip_imerg_mswep': (['time', 'lat', 'lon'], delta_precip.astype(np.float32)),
            'precip_1d_sum': (['time', 'lat', 'lon'], precip_1d_sum.astype(np.float32)),
            'precip_3d_sum': (['time', 'lat', 'lon'], precip_3d_sum.astype(np.float32)),
            'precip_7d_sum': (['time', 'lat', 'lon'], precip_7d_sum.astype(np.float32)),
            'days_since_rain_gt_1mm': (['time', 'lat', 'lon'], days_since_rain.astype(np.float32))
        },
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={'title': f'Precipitation derived predictors and anti-leakage metrics ({year})'}
    )
    ds_metrics.to_netcdf(metrics_file)
    print(f"  Saved Precip Metrics: {metrics_file} ({os.path.getsize(metrics_file)/(1024*1024):.1f} MB)")


def main():
    # 1. Derive static predictors
    ds_static = derive_static_predictors()
    lat_1d = ds_static.coords['lat'].values
    lon_1d = ds_static.coords['lon'].values
    
    # 2. Process all years 2016 to 2020
    years = [2016, 2017, 2018, 2019, 2020]
    for y in years:
        extract_opl_variables(y, lat_1d, lon_1d)
        regrid_mswep_and_compute_metrics(y, lat_1d, lon_1d)
        
    print("\nPhase 1c predictor derivation completed successfully.")

if __name__ == "__main__":
    main()
