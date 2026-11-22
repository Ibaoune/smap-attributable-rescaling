"""
SAR I/O utilities: reading, writing, configuration parsing, and provenance tracking.
"""

import os
import hashlib
import datetime
import yaml
import numpy as np
import pandas as pd
import xarray as xr


def compute_sha256(filepath: str, chunk_size: int = 65536) -> str:
    """Compute the SHA256 checksum of a file."""
    sha = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            sha.update(data)
    return sha.hexdigest()


def log_provenance(manifest_csv: str, source_path: str, target_path: str,
                   producing_script: str, description: str = "") -> None:
    """
    Log an entry into the provenance manifest CSV.
    """
    size_bytes = os.path.getsize(source_path)
    sha256 = compute_sha256(source_path)
    mtime = os.path.getmtime(source_path)
    source_date = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')

    record = {
        'source_path': source_path,
        'target_path': target_path,
        'size_bytes': size_bytes,
        'sha256': sha256,
        'source_date': source_date,
        'producing_script': producing_script,
        'description': description
    }

    df_record = pd.DataFrame([record])
    if os.path.exists(manifest_csv) and os.path.getsize(manifest_csv) > 0:
        df_record.to_csv(manifest_csv, mode='a', header=False, index=False)
    else:
        df_record.to_csv(manifest_csv, mode='w', header=True, index=False)


def load_yaml(filepath: str) -> dict:
    """Load and parse a YAML configuration file."""
    with open(filepath, 'r') as f:
        return yaml.safe_load(f)


def load_canonical_grid(ref_nc_path: str):
    """
    Load canonical NorthMor 121 x 181 coordinates and land mask.
    Returns (lat_1d, lon_1d, land_mask_2d).
    """
    ds = xr.open_dataset(ref_nc_path)
    lat = ds['lat'].values
    lon = ds['lon'].values
    
    # Identify land mask from observation_frequency or ssm
    if 'observation_frequency' in ds:
        mask = ~np.isnan(ds['observation_frequency'].values[0])
    elif 'ssm' in ds:
        mask = ~np.isnan(ds['ssm'].values[0])
    else:
        first_var = list(ds.data_vars.keys())[0]
        val = ds[first_var].values
        mask = ~np.isnan(val[0] if val.ndim == 3 else val)
    ds.close()
    return lat, lon, mask
