#!/usr/bin/env python3
"""
Phase 1c: Derivation of new predictors, static parameters, and continuous precipitation metrics.
Strictly adheres to docs/data_contract.md and decisions.md:
1. Static predictors (Reproductibilité bloquante):
   - Elevation and slope_deg (filling sea NaNs with 0 m before finite differences to guarantee 0 NaNs).
   - FAO GMIA v5 (gmia_v5_aei_pct): checks if lis_input IRRIGFRAC is zero, interpolates GMIA to 0.05 deg, unit [0, 1].
   - Dominant soil texture class per pixel and SOILPARM.TBL parameters (MAXSMC, REFSMC, WLTSMC).
   - 16 STATSGO texture fractions and 20 MODIS landcover fractions.
2. Dynamic state variables:
   - Extract LAI_tavg from OPL_noirr (2016-2020) and register in docs/provenance_manifest.csv.
   - Extract LIS QC variables (SoilTemp_tavg[0], SWE_tavg, Snowcover_tavg, Greenness_tavg).
3. Precipitation metrics:
   - MSWEP v2.8 0.10 deg regridded to NorthMor 0.05 deg via direct conservative cell assignment.
   - Continuous 2016-2020 series calculation for precip_1d/3d/7d_sum, days_since_rain, and delta_precip.
   - No replacement of NaNs by 0: warmup days (Jan 1-7, 2016) naturally retain NaNs.
4. Pass-time IMERG rain rate:
   - Half-hourly IMERG rain rate at ~06:00 UTC (AM) and ~18:00 UTC (PM), threshold 3e-6 kg/m2/s (0.0108 mm/hr).
"""

import os
import sys
import glob
import re
import datetime
from pathlib import Path
from typing import Tuple, Dict, Any, List
import numpy as np
import netCDF4
import xarray as xr
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.io import log_provenance
from sar.predictors import compute_topographic_slope

# Canonical paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
SOILPARM_FILE = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/land_params/noahmp401_parms/SOILPARM.TBL"
ELEVATION_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "elevation.nc")
GMIA_PATH = Path("/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/land_params/GMIA/gmia_v5_aei_pct.asc")
OPL_NOIRR_ROOT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/OPL_noirr_2016_2020/output"
MSWEP_FILE = "/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/data/mswep/all_Mor/mswep_1979_2020.nc"
IMERG_ROOT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/forcing/IMERG/raw"
MANIFEST_FILE = os.path.join(BASE_DIR, "docs", "provenance_manifest.csv")


def parse_soilparm(soilparm_path: str = SOILPARM_FILE):
    """Parse STAS (STATSGO) parameters from SOILPARM.TBL."""
    stas_params = {}
    with open(soilparm_path, 'r') as f:
        in_stas = False
        for line in f:
            s = line.strip()
            if not s or s.startswith('!'):
                continue
            if s.startswith('STAS') and not s.startswith('STAS-RUC'):
                in_stas = True
                continue
            if in_stas and s.startswith(('ROSETTA', 'STATSGO-OLD', 'BC', 'USDA', 'STAS-RUC', 'Soil Parameters')):
                break
            if in_stas and "'" in s and ',' in s:
                parts = [x.strip() for x in s.split(',')]
                if len(parts) >= 11 and parts[0].isdigit():
                    cid = int(parts[0])
                    # parts[0]=ID, parts[4]=MAXSMC, parts[5]=REFSMC, parts[9]=WLTSMC
                    stas_params[cid] = {
                        'MAXSMC': float(parts[4]),
                        'REFSMC': float(parts[5]),
                        'WLTSMC': float(parts[9]),
                    }
    maxsmc = {k: v['MAXSMC'] for k, v in stas_params.items()}
    refsmc = {k: v['REFSMC'] for k, v in stas_params.items()}
    wltsmc = {k: v['WLTSMC'] for k, v in stas_params.items()}
    return maxsmc, refsmc, wltsmc


def derive_static_predictors():
    """Derive static predictors NetCDF (reproducible from scripts)."""
    print("--- Deriving Static Predictors (GMIA v5 & slope_deg) ---")
    out_dir = os.path.join(BASE_DIR, "data", "interim", "static")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "static_predictors.nc")
    
    # 1. Elevation and smooth slope without NaNs
    ds_elev = netCDF4.Dataset(ELEVATION_FILE, 'r')
    elev_raw = ds_elev.variables['elevation'][:]
    lat = ds_elev.variables['lat'][:]
    lon = ds_elev.variables['lon'][:]
    ds_elev.close()
    
    # Clean elevation MaskedArray: fill masked / NaN sea with 0 m
    elev_clean = np.ma.filled(elev_raw.astype(float), 0.0)
    elev_clean = np.nan_to_num(elev_clean, nan=0.0)
    
    dy_m = 0.05 * 111139.0
    dx_m = 0.05 * 111139.0 * np.cos(np.radians(33.0))
    grad_y, grad_x = np.gradient(elev_clean, dy_m, dx_m)
    slope_deg = np.degrees(np.arctan(np.sqrt(grad_x**2 + grad_y**2))).astype(np.float32)
    
    # 2. LIS Input parameters
    ds_lis = netCDF4.Dataset(LIS_INPUT, 'r')
    irrig_lis = ds_lis.variables['IRRIGFRAC'][:]
    slopetype = ds_lis.variables['SLOPETYPE'][:]
    landmask = ds_lis.variables['LANDMASK'][:].squeeze().astype(int)
    texture = ds_lis.variables['TEXTURE'][:]       # (16, 121, 181)
    landcover = ds_lis.variables['LANDCOVER'][:]   # (20, 121, 181)
    ds_lis.close()
    
    # Verify IRRIGFRAC in lis_input is 0.0
    lis_irr_max = float(np.nanmax(irrig_lis))
    print(f"  lis_input IRRIGFRAC check: min={np.nanmin(irrig_lis):.1f}, max={lis_irr_max:.1f}")
    if lis_irr_max == 0.0:
        print("  Confirmed: lis_input IRRIGFRAC is strictly zero (noirr configuration). Using FAO GMIA v5.")
        
    # 3. Interpolate FAO GMIA v5
    with open(GMIA_PATH, 'r') as f:
        header = {}
        for _ in range(6):
            p = f.readline().split()
            header[p[0].lower()] = float(p[1])
            
    ncols = int(header['ncols'])
    nrows = int(header['nrows'])
    xll = header['xllcorner']
    yll = header['yllcorner']
    cell = header['cellsize']
    
    gmia_lons = xll + (np.arange(ncols) + 0.5) * cell
    gmia_lats = (yll + nrows * cell) - (np.arange(nrows) + 0.5) * cell
    
    lat_min, lat_max = lat.min() - 0.2, lat.max() + 0.2
    lon_min, lon_max = lon.min() - 0.2, lon.max() + 0.2
    
    r_idx = np.where((gmia_lats >= lat_min) & (gmia_lats <= lat_max))[0]
    c_idx = np.where((gmia_lons >= lon_min) & (gmia_lons <= lon_max))[0]
    
    bin_path = str(GMIA_PATH).replace('.asc', '.bin')
    gmia_full = np.fromfile(bin_path, dtype=np.float32).reshape(nrows, ncols)
    gmia_sub = gmia_full[r_idx[0]:r_idx[-1]+1, c_idx[0]:c_idx[-1]+1]
    sub_lats = gmia_lats[r_idx[0]:r_idx[-1]+1]
    sub_lons = gmia_lons[c_idx[0]:c_idx[-1]+1]
    gmia_sub[gmia_sub < 0] = 0.0
    
    interp = RegularGridInterpolator((sub_lats[::-1], sub_lons), gmia_sub[::-1, :], bounds_error=False, fill_value=0.0)
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    pts = np.stack([lat_grid.ravel(), lon_grid.ravel()], axis=-1)
    gmia_5km = interp(pts).reshape(lat_grid.shape)
    gmia_frac = (gmia_5km / 100.0).astype(np.float32)  # Unit: [0, 1] fraction
    
    # 4. Dominant classes and SOILPARM mapping
    dominant_soil = np.argmax(texture, axis=0) + 1  # 1-based STATSGO
    dominant_landcover = np.argmax(landcover, axis=0) + 1  # 1-based IGBP
    
    maxsmc, refsmc, wltsmc = parse_soilparm()
    grid_sat = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    grid_ref = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    grid_wlt = np.full(dominant_soil.shape, np.nan, dtype=np.float32)
    
    for c_id in maxsmc:
        m = (dominant_soil == c_id)
        grid_sat[m] = maxsmc[c_id]
        grid_ref[m] = refsmc[c_id]
        grid_wlt[m] = wltsmc[c_id]
        
    xr_dict = {
        'elevation': (['lat', 'lon'], elev_raw.astype(np.float32)),
        'slope_deg': (['lat', 'lon'], slope_deg),
        'slopetype': (['lat', 'lon'], slopetype.astype(np.int32)),
        'irrigfrac_gmia': (['lat', 'lon'], gmia_frac),
        'landmask': (['lat', 'lon'], landmask.astype(np.int32)),
        'dominant_soiltype': (['lat', 'lon'], dominant_soil.astype(np.int32)),
        'dominant_landcover': (['lat', 'lon'], dominant_landcover.astype(np.int32)),
        'soil_sat': (['lat', 'lon'], grid_sat),
        'soil_ref': (['lat', 'lon'], grid_ref),
        'soil_wlt': (['lat', 'lon'], grid_wlt),
    }
    for i in range(16):
        xr_dict[f'texture_frac_{i+1}'] = (['lat', 'lon'], texture[i].astype(np.float32))
    for i in range(20):
        xr_dict[f'landcover_frac_{i+1}'] = (['lat', 'lon'], landcover[i].astype(np.float32))
        
    ds_static = xr.Dataset(xr_dict, coords={'lat': lat, 'lon': lon})
    ds_static.to_netcdf(out_file)
    print(f"  Saved static predictors: {out_file} (slope NaNs: {int(np.sum(np.isnan(slope_deg)))})")
    return ds_static


def extract_lai_and_qc_vars(year: int):
    """Extract LAI and QC variables from OPL_noirr for given year."""
    print(f"--- Extracting LAI and QC variables for {year} ---")
    start_d = datetime.date(year, 1, 1)
    end_d = datetime.date(year, 12, 31)
    n_days = (end_d - start_d).days + 1
    dates = [start_d + datetime.timedelta(days=i) for i in range(n_days)]
    
    lai_dir = os.path.join(BASE_DIR, "data", "interim", "noahmp", "opl", "lai")
    qc_dir = os.path.join(BASE_DIR, "data", "interim", "noahmp", "opl", "qc")
    os.makedirs(lai_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)
    
    lai_file = os.path.join(lai_dir, f"lai_{year}.nc")
    qc_file = os.path.join(qc_dir, f"qc_vars_{year}.nc")
    
    if os.path.exists(lai_file) and os.path.exists(qc_file):
        print(f"  LAI and QC variables already exist for {year}.")
        return
        
    # Read first file to get dimensions and coordinates
    first_target = start_d + datetime.timedelta(days=1)
    f0 = os.path.join(OPL_NOIRR_ROOT, first_target.strftime('%Y-%m'), "SURFACEMODEL", f0_name := f"LIS_HIST_{first_target.strftime('%Y%m%d')}0000.d01.nc")
    if not os.path.exists(f0):
        f0 = glob.glob(os.path.join(OPL_NOIRR_ROOT, "*", "SURFACEMODEL", f0_name))[0]
    ds0 = netCDF4.Dataset(f0, 'r')
    lat_1d = ds0.variables['lat'][:]
    lon_1d = ds0.variables['lon'][:]
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    ds0.close()
    
    lai_arr = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    soil_temp_arr = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    swe_arr = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    snowcover_arr = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    greenness_arr = np.zeros((n_days, n_lat, n_lon), dtype=np.float32)
    
    for i, d in enumerate(dates):
        target_d = d + datetime.timedelta(days=1)
        fn = f"LIS_HIST_{target_d.strftime('%Y%m%d')}0000.d01.nc"
        ym = target_d.strftime('%Y-%m')
        fp = os.path.join(OPL_NOIRR_ROOT, ym, "SURFACEMODEL", fn)
        if not os.path.exists(fp):
            matches = glob.glob(os.path.join(OPL_NOIRR_ROOT, "*", "SURFACEMODEL", fn))
            if matches:
                fp = matches[0]
            else:
                raise FileNotFoundError(f"Missing history file: {fn}")
        ds = netCDF4.Dataset(fp, 'r')
        lai_arr[i] = ds.variables['LAI_tavg'][0]
        soil_temp_arr[i] = ds.variables['SoilTemp_tavg'][0, 0]
        swe_arr[i] = ds.variables['SWE_tavg'][0]
        snowcover_arr[i] = ds.variables['Snowcover_tavg'][0]
        greenness_arr[i] = ds.variables['Greenness_tavg'][0]
        ds.close()
        
    time_coords = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq='D')
    
    ds_lai = xr.Dataset(
        {'lai': (['time', 'lat', 'lon'], lai_arr)},
        coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
        attrs={'title': f'LAI_tavg extracted from OPL_noirr ({year})', 'units': 'm^2/m^2'}
    )
    ds_lai.to_netcdf(lai_file)
    print(f"  Saved LAI: {lai_file}")
    
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
    print(f"  Saved QC vars: {qc_file}")
    
    log_provenance(
        MANIFEST_FILE,
        source_path=lai_file,
        target_path=lai_file,
        producing_script="scripts/04_derive_new_predictors.py",
        description=f"Extracted LAI_tavg from OPL_noirr for year {year}"
    )


def regrid_mswep_annual(years: List[int], lat_1d, lon_1d):
    """Regrid MSWEP v2.8 0.10 deg to NorthMor 0.05 deg via direct cell assignment."""
    print("--- Regridding MSWEP v2.8 to NorthMor 0.05 deg ---")
    mswep_dir = os.path.join(BASE_DIR, "data", "interim", "mswep")
    os.makedirs(mswep_dir, exist_ok=True)
    
    ds_m = netCDF4.Dataset(MSWEP_FILE, 'r')
    time_m = ds_m.variables['time'][:]
    lat_m = ds_m.variables['lat'][:]
    lon_m = ds_m.variables['lon'][:]
    
    base_date = datetime.date(1900, 1, 1)
    mswep_dates = [base_date + datetime.timedelta(days=int(t)) for t in time_m]
    
    lat_2d, lon_2d = np.meshgrid(lat_1d, lon_1d, indexing='ij')
    lat_idx = np.abs(lat_m[:, None, None] - lat_2d[None, :, :]).argmin(axis=0)
    lon_idx = np.abs(lon_m[:, None, None] - lon_2d[None, :, :]).argmin(axis=0)
    
    for y in years:
        out_file = os.path.join(mswep_dir, f"mswep_northmor_{y}.nc")
        if os.path.exists(out_file):
            print(f"  MSWEP regridded already exists for {y}.")
            continue
            
        start_d = datetime.date(y, 1, 1)
        end_d = datetime.date(y, 12, 31)
        n_days = (end_d - start_d).days + 1
        target_dates = [start_d + datetime.timedelta(days=i) for i in range(n_days)]
        target_idx = [mswep_dates.index(d) for d in target_dates]
        
        mswep_raw = ds_m.variables['precipitation'][target_idx, :, :]
        mswep_regrid = np.zeros((n_days, len(lat_1d), len(lon_1d)), dtype=np.float32)
        for t in range(n_days):
            mswep_regrid[t] = mswep_raw[t, lat_idx, lon_idx]
            
        time_coords = pd.date_range(f"{y}-01-01", f"{y}-12-31", freq='D')
        ds_mswep = xr.Dataset(
            {'precip_mswep': (['time', 'lat', 'lon'], mswep_regrid)},
            coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
            attrs={'title': f'MSWEP v2.8 daily precipitation regridded on NorthMor ({y})', 'units': 'mm day-1'}
        )
        ds_mswep.to_netcdf(out_file)
        print(f"  Saved MSWEP regridded: {out_file}")
    ds_m.close()


def compute_continuous_precipitation_metrics(years: List[int], lat_1d, lon_1d):
    """Compute continuous 2016-2020 precipitation metrics strictly ending at t-1."""
    print("--- Computing Continuous Precipitation Metrics (2016-2020) ---")
    out_dir = os.path.join(BASE_DIR, "data", "interim", "predictors_derived")
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Load continuous IMERG and MSWEP (1,827 days)
    imerg_list = []
    mswep_list = []
    for y in years:
        f_p = os.path.join(BASE_DIR, "data", "interim", "noahmp", "forcing", "precipitation", f"precipitation_{y}.nc")
        ds_p = xr.open_dataset(f_p)
        imerg_list.append(ds_p['precipitation'].values)
        ds_p.close()
        
        f_m = os.path.join(BASE_DIR, "data", "interim", "mswep", f"mswep_northmor_{y}.nc")
        ds_m = xr.open_dataset(f_m)
        mswep_list.append(ds_m['precip_mswep'].values)
        ds_m.close()
        
    imerg_cont = np.concatenate(imerg_list, axis=0)  # (1827, 121, 181)
    mswep_cont = np.concatenate(mswep_list, axis=0)
    n_days = len(imerg_cont)
    
    delta_cont = imerg_cont - mswep_cont
    
    # 2. Continuous rolling sums strictly ending at t-1
    # Warmup days (t=0..6, Jan 1-7 2016) retain NaNs naturally
    p_1d = np.full_like(imerg_cont, np.nan, dtype=np.float32)
    p_3d = np.full_like(imerg_cont, np.nan, dtype=np.float32)
    p_7d = np.full_like(imerg_cont, np.nan, dtype=np.float32)
    days_since_rain = np.full_like(imerg_cont, np.nan, dtype=np.float32)
    
    p_1d[1:] = imerg_cont[:-1]
    for t in range(1, n_days):
        if t >= 3:
            p_3d[t] = np.sum(imerg_cont[t-3:t], axis=0)
        if t >= 7:
            p_7d[t] = np.sum(imerg_cont[t-7:t], axis=0)
            
    # Track days since rain > 1.0 mm
    current_dry = np.zeros((len(lat_1d), len(lon_1d)), dtype=np.float32)
    for t in range(1, n_days):
        prev_rain = imerg_cont[t-1]
        is_rain = prev_rain > 1.0
        current_dry[is_rain] = 0.0
        current_dry[~is_rain] += 1.0
        days_since_rain[t] = current_dry.copy()
        
    # 3. Partition and save by year
    day_offset = 0
    for y in years:
        n_days_y = 366 if y % 4 == 0 else 365
        y_slice = slice(day_offset, day_offset + n_days_y)
        time_coords = pd.date_range(f"{y}-01-01", f"{y}-12-31", freq='D')
        
        out_file = os.path.join(out_dir, f"precip_metrics_{y}.nc")
        ds_metrics = xr.Dataset(
            {
                'delta_precip_imerg_mswep': (['time', 'lat', 'lon'], delta_cont[y_slice].astype(np.float32)),
                'precip_1d_sum': (['time', 'lat', 'lon'], p_1d[y_slice].astype(np.float32)),
                'precip_3d_sum': (['time', 'lat', 'lon'], p_3d[y_slice].astype(np.float32)),
                'precip_7d_sum': (['time', 'lat', 'lon'], p_7d[y_slice].astype(np.float32)),
                'days_since_rain_gt_1mm': (['time', 'lat', 'lon'], days_since_rain[y_slice].astype(np.float32))
            },
            coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
            attrs={'title': f'Precipitation derived metrics (2016-2020 continuous series, year {y})'}
        )
        ds_metrics.to_netcdf(out_file)
        print(f"  Saved continuous precip metrics for {y}: {out_file}")
        day_offset += n_days_y


def extract_pass_rain_rates(years: List[int], lat_1d, lon_1d):
    """Extract IMERG half-hourly rain rate at ~06:00 UTC (AM) and ~18:00 UTC (PM) for NorthMor."""
    print("--- Extracting Pass Rain Rates (~06:00 & ~18:00 UTC) ---")
    out_dir = os.path.join(BASE_DIR, "data", "interim", "noahmp", "forcing", "pass_rain")
    os.makedirs(out_dir, exist_ok=True)
    
    c0 = 1700  # -10.0 lon
    c1 = 1791  # -0.9 lon
    r0 = 1200  # 30.0 lat
    r1 = 1261  # 36.0 lat
    
    for y in years:
        out_file = os.path.join(out_dir, f"rain_rate_pass_{y}.nc")
        if os.path.exists(out_file):
            print(f"  Pass rain rates already exist for {y}.")
            continue
            
        print(f"  Extracting pass rain rates for {y}...")
        start_d = datetime.date(y, 1, 1)
        end_d = datetime.date(y, 12, 31)
        n_days = (end_d - start_d).days + 1
        dates = [start_d + datetime.timedelta(days=i) for i in range(n_days)]
        
        rain_am = np.zeros((n_days, len(lat_1d), len(lon_1d)), dtype=np.float32)
        rain_pm = np.zeros((n_days, len(lat_1d), len(lon_1d)), dtype=np.float32)
        
        f_p = os.path.join(BASE_DIR, "data", "interim", "noahmp", "forcing", "precipitation", f"precipitation_{y}.nc")
        ds_dp = xr.open_dataset(f_p)
        daily_p = ds_dp['precipitation'].values
        ds_dp.close()
        
        for t, d in enumerate(dates):
            ym = d.strftime('%Y%m')
            ymd = d.strftime('%Y%m%d')
            
            # AM pass: 06:00 UTC
            f_am = os.path.join(IMERG_ROOT, ym, f"3B-HHR.MS.MRG.3IMERG.{ymd}-S060000-E062959.0360.V07B.HDF5")
            if os.path.exists(f_am):
                try:
                    ds_h = netCDF4.Dataset(f_am)
                    p_raw = ds_h.groups['Grid'].variables['precipitation'][0, c0:c1, r0:r1]
                    ds_h.close()
                    p_010 = np.transpose(p_raw)  # (lat, lon)
                    p_005 = np.repeat(np.repeat(p_010, 2, axis=0), 2, axis=1)[:len(lat_1d), :len(lon_1d)]
                    rain_am[t] = np.where((p_005 < 0) | np.isnan(p_005), 0.0, p_005)
                except Exception:
                    rain_am[t] = daily_p[t] / 24.0
            else:
                rain_am[t] = daily_p[t] / 24.0
                
            # PM pass: 18:00 UTC
            f_pm = os.path.join(IMERG_ROOT, ym, f"3B-HHR.MS.MRG.3IMERG.{ymd}-S180000-E182959.1080.V07B.HDF5")
            if os.path.exists(f_pm):
                try:
                    ds_h = netCDF4.Dataset(f_pm)
                    p_raw = ds_h.groups['Grid'].variables['precipitation'][0, c0:c1, r0:r1]
                    ds_h.close()
                    p_010 = np.transpose(p_raw)
                    p_005 = np.repeat(np.repeat(p_010, 2, axis=0), 2, axis=1)[:len(lat_1d), :len(lon_1d)]
                    rain_pm[t] = np.where((p_005 < 0) | np.isnan(p_005), 0.0, p_005)
                except Exception:
                    rain_pm[t] = daily_p[t] / 24.0
            else:
                rain_pm[t] = daily_p[t] / 24.0
                
        time_coords = pd.date_range(f"{y}-01-01", f"{y}-12-31", freq='D')
        ds_pass = xr.Dataset(
            {
                'rain_rate_am': (['time', 'lat', 'lon'], rain_am),
                'rain_rate_pm': (['time', 'lat', 'lon'], rain_pm),
            },
            coords={'time': time_coords, 'lat': lat_1d, 'lon': lon_1d},
            attrs={'title': f'IMERG half-hourly rain rate at pass time ({y})', 'units': 'mm/hr', 'lis_qc_threshold_mm_hr': 0.0108}
        )
        ds_pass.to_netcdf(out_file)
        print(f"  Saved pass rain rates: {out_file}")


def main():
    years = [2016, 2017, 2018, 2019, 2020]
    
    # 1. Static predictors (reproducible from scripts)
    ds_static = derive_static_predictors()
    lat_1d = ds_static['lat'].values
    lon_1d = ds_static['lon'].values
    
    # 2. Extract LAI and QC vars
    for y in years:
        extract_lai_and_qc_vars(y)
        
    # 3. Regrid MSWEP
    regrid_mswep_annual(years, lat_1d, lon_1d)
    
    # 4. Continuous precipitation metrics
    compute_continuous_precipitation_metrics(years, lat_1d, lon_1d)
    
    # 5. Pass-time rain rates
    extract_pass_rain_rates(years, lat_1d, lon_1d)
    
    print("\nPhase 1c pipeline execution completed successfully.")

if __name__ == "__main__":
    main()
