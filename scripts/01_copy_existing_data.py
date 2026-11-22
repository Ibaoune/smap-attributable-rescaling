#!/usr/bin/env python3
"""
Step 1a: Copy and verify existing data from smap-xai-audit to smap-attributable-rescaling.
Registers provenance manifest and validates spatial and temporal integrity.
"""

import os
import shutil
import hashlib
import datetime
import numpy as np
import xarray as xr
import pandas as pd

SRC_ROOT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-xai-audit"
DEST_ROOT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling"
MANIFEST_FILE = os.path.join(DEST_ROOT, "docs", "provenance_manifest.csv")

YEARS = [2016, 2017, 2018, 2019, 2020]
EXPECTED_DAYS = {2016: 366, 2017: 365, 2018: 365, 2019: 365, 2020: 366}

def compute_sha256(filepath):
    sha = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(65536)
            if not data:
                break
            sha.update(data)
    return sha.hexdigest()

def main():
    print("=== Step 1a: Starting copy and integrity verification ===")
    
    files_to_copy = []

    # 1. Statiques
    statics = [
        ("data/processed/static/elevation/elevation.nc", "21_prepare_elevation.py", "Static surface elevation"),
        ("data/processed/static/soil_texture/soil_texture_fraction.nc", "22_prepare_soil_texture.py", "Static 16 STATSGO soil texture fractions"),
        ("data/processed/static/land_cover/land_cover_fraction.nc", "23_prepare_land_cover.py", "Static 20 MODIS-IGBP land cover fractions")
    ]
    for rel_path, script, desc in statics:
        src = os.path.join(SRC_ROOT, rel_path)
        dst = os.path.join(DEST_ROOT, "data", "raw_copied", os.path.basename(rel_path))
        files_to_copy.append({
            'src': src, 'dst': dst, 'script': script, 'desc': desc, 'is_static': True
        })

    # 2. OL variables (2016-2020)
    ol_vars = [
        ("ssm", "data/interim/noahmp/opl/ssm/ssm_{year}.nc", "01_prepare_ssm.py", "Noah-MP OL surface soil moisture layer 1 daily"),
        ("sm_layer2", "data/interim/noahmp/opl/sm_layer2/sm_layer2_{year}.nc", "02_prepare_soil_moisture_layer.py", "Noah-MP OL soil moisture layer 2 daily"),
        ("sm_layer3", "data/interim/noahmp/opl/sm_layer3/sm_layer3_{year}.nc", "02_prepare_soil_moisture_layer.py", "Noah-MP OL soil moisture layer 3 daily"),
        ("sm_layer4", "data/interim/noahmp/opl/sm_layer4/sm_layer4_{year}.nc", "02_prepare_soil_moisture_layer.py", "Noah-MP OL soil moisture layer 4 daily"),
        ("rzsm", "data/processed/daily_model/opl/rzsm/rzsm_{year}.nc", "03_derive_rzsm.py", "Noah-MP OL root zone soil moisture 0-100cm daily"),
        ("surface_runoff", "data/interim/noahmp/opl/surface_runoff/surface_runoff_{year}.nc", "07_prepare_surface_runoff.py", "Noah-MP OL surface runoff daily"),
        ("baseflow", "data/interim/noahmp/opl/baseflow/baseflow_{year}.nc", "08_prepare_baseflow.py", "Noah-MP OL baseflow daily"),
        ("et", "data/interim/noahmp/opl/et/et_{year}.nc", "06_prepare_et.py", "Noah-MP OL evapotranspiration daily"),
    ]
    for var, template, script, desc in ol_vars:
        for y in YEARS:
            rel_src = template.format(year=y)
            src = os.path.join(SRC_ROOT, rel_src)
            dst_sub = os.path.join(DEST_ROOT, "data", "interim", "noahmp", "opl", var)
            os.makedirs(dst_sub, exist_ok=True)
            dst = os.path.join(dst_sub, f"{var}_{y}.nc")
            files_to_copy.append({
                'src': src, 'dst': dst, 'script': script, 'desc': f"{desc} ({y})", 'is_static': False, 'year': y, 'var': var
            })

    # 3. Forçages (2016-2020)
    forcing_vars = [
        ("precipitation", "data/interim/noahmp/forcing/precipitation/precipitation_{year}.nc", "04_prepare_precipitation.py", "GPM IMERG V07B precipitation daily"),
        ("air_temperature", "data/interim/noahmp/forcing/air_temperature/air_temperature_{year}.nc", "05_prepare_air_temperature.py", "MERRA-2 2m air temperature daily mean"),
    ]
    for var, template, script, desc in forcing_vars:
        for y in YEARS:
            rel_src = template.format(year=y)
            src = os.path.join(SRC_ROOT, rel_src)
            dst_sub = os.path.join(DEST_ROOT, "data", "interim", "noahmp", "forcing", var)
            os.makedirs(dst_sub, exist_ok=True)
            dst = os.path.join(dst_sub, f"{var}_{y}.nc")
            files_to_copy.append({
                'src': src, 'dst': dst, 'script': script, 'desc': f"{desc} ({y})", 'is_static': False, 'year': y, 'var': var
            })

    # 4. Diagnostics EnKF (2016-2020)
    # da_nocdf
    enkf_nocdf = [
        ("daily_innovation", "data/processed/daily_enkf/da_nocdf/daily_innovation/daily_innovation_{year}.nc", "17_derive_daily_innovation.py", "Daily innovation da_nocdf"),
        ("daily_signed_increment", "data/processed/daily_enkf/da_nocdf/daily_signed_increment/daily_signed_increment_{year}.nc", "15_derive_daily_signed_increment.py", "Daily signed increment da_nocdf")
    ]
    for var, template, script, desc in enkf_nocdf:
        for y in YEARS:
            rel_src = template.format(year=y)
            src = os.path.join(SRC_ROOT, rel_src)
            dst_sub = os.path.join(DEST_ROOT, "data", "processed", "daily_enkf", "da_nocdf", var)
            os.makedirs(dst_sub, exist_ok=True)
            dst = os.path.join(dst_sub, f"{var}_{y}.nc")
            files_to_copy.append({
                'src': src, 'dst': dst, 'script': script, 'desc': f"{desc} ({y})", 'is_static': False, 'year': y, 'var': var, 'check_cdf_date': False
            })

    # da_cdf (CDF v2 uniquement, verifie > 03/11/2026)
    enkf_cdf = [
        ("daily_innovation", "data/processed/daily_enkf/da_cdf/daily_innovation/daily_innovation_{year}.nc", "17_derive_daily_innovation.py", "Daily innovation da_cdf v2"),
        ("daily_signed_increment", "data/processed/daily_enkf/da_cdf/daily_signed_increment/daily_signed_increment_{year}.nc", "15_derive_daily_signed_increment.py", "Daily signed increment da_cdf v2")
    ]
    for var, template, script, desc in enkf_cdf:
        for y in YEARS:
            rel_src = template.format(year=y)
            src = os.path.join(SRC_ROOT, rel_src)
            dst_sub = os.path.join(DEST_ROOT, "data", "processed", "daily_enkf", "da_cdf", var)
            os.makedirs(dst_sub, exist_ok=True)
            dst = os.path.join(dst_sub, f"{var}_{y}.nc")
            files_to_copy.append({
                'src': src, 'dst': dst, 'script': script, 'desc': f"{desc} ({y})", 'is_static': False, 'year': y, 'var': var, 'check_cdf_date': True
            })

    print(f"Total files identified to copy: {len(files_to_copy)}")

    # CDF v2 threshold: 2026-11-03 00:00:00
    cdf_v2_cutoff = datetime.datetime(2026, 11, 3, 0, 0, 0).timestamp()

    manifest_rows = []
    
    # Execute copy, hashing, and validation
    for i, item in enumerate(files_to_copy, 1):
        src = item['src']
        dst = item['dst']
        if not os.path.exists(src):
            raise FileNotFoundError(f"Source file missing: {src}")

        # Check date for da_cdf
        mtime = os.path.getmtime(src)
        mtime_dt = datetime.datetime.fromtimestamp(mtime)
        if item.get('check_cdf_date', False):
            if mtime < cdf_v2_cutoff:
                raise ValueError(f"CDF v2 validation failed: {src} has mtime {mtime_dt} which is older than 03/11/2026!")

        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)

        size_bytes = os.path.getsize(dst)
        sha256 = compute_sha256(dst)
        source_date = mtime_dt.strftime('%Y-%m-%d %H:%M:%S')

        manifest_rows.append({
            'source_path': src,
            'target_path': dst,
            'size_bytes': size_bytes,
            'sha256': sha256,
            'source_date': source_date,
            'producing_script': item['script'],
            'description': item['desc']
        })

        # Integrity check via xarray
        with xr.open_dataset(dst) as ds:
            # Lat/Lon dimensions
            if 'lat' in ds.coords and 'lon' in ds.coords:
                n_lat = len(ds['lat'])
                n_lon = len(ds['lon'])
            elif 'lat' in ds.dims and 'lon' in ds.dims:
                n_lat = ds.dims['lat']
                n_lon = ds.dims['lon']
            else:
                n_lat = ds[list(ds.data_vars)[0]].shape[-2]
                n_lon = ds[list(ds.data_vars)[0]].shape[-1]
                
            if n_lat != 121 or n_lon != 181:
                raise ValueError(f"Integrity check failed for {dst}: shape is ({n_lat}, {n_lon}), expected (121, 181)")

            # Time dimension check
            if not item['is_static']:
                n_time = len(ds['time'])
                expected = EXPECTED_DAYS[item['year']]
                if n_time != expected:
                    raise ValueError(f"Integrity check failed for {dst}: {n_time} days found, expected {expected}")

        if i % 10 == 0 or i == len(files_to_copy):
            print(f"[{i:02d}/{len(files_to_copy)}] Copied and verified: {os.path.basename(dst)}")

    # Write manifest
    df_manifest = pd.DataFrame(manifest_rows)
    df_manifest.to_csv(MANIFEST_FILE, index=False)
    print(f"\nManifest successfully updated with {len(df_manifest)} entries in {MANIFEST_FILE}")
    print("=== Step 1a completed successfully! ===")

if __name__ == "__main__":
    main()
