"""
SMAP SPL3SMP_E Raw HDF5 Data Reader.
Handles extraction of AM and PM retrieval fields, quality flags, and geophysical variables.
Standardizes PM variable names (which carry the '_pm' suffix in HDF5) to canonical variable names.
"""

import os
from typing import Dict, Optional, Tuple
import numpy as np
import netCDF4


SMAP_CANONICAL_VARIABLES = [
    'soil_moisture',
    'retrieval_qual_flag',
    'surface_flag',
    'vegetation_opacity',
    'vegetation_water_content',
    'surface_temperature',
    'roughness_coefficient',
    'static_water_body_fraction',
]

GROUP_AM = 'Soil_Moisture_Retrieval_Data_AM'
GROUP_PM = 'Soil_Moisture_Retrieval_Data_PM'


def open_smap_h5(filepath: str) -> netCDF4.Dataset:
    """Open a raw SMAP SPL3SMP_E HDF5 file using netCDF4."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"SMAP raw file not found: {filepath}")
    return netCDF4.Dataset(filepath, 'r')


def read_smap_pass_data(ds: netCDF4.Dataset, pass_name: str = 'AM') -> Dict[str, np.ndarray]:
    """
    Extract all required SMAP variables for a specific overpass (AM or PM).
    Standardizes variable names (removes '_pm' suffix for PM pass).
    Returns a dictionary of 2D numpy arrays on the global EASE-Grid 2.0 (1624 x 3856).
    """
    pass_upper = pass_name.upper()
    group_name = GROUP_AM if pass_upper == 'AM' else GROUP_PM
    if group_name not in ds.groups:
        raise KeyError(f"Group {group_name} not found in SMAP file.")

    grp = ds.groups[group_name]
    data = {}
    for var in SMAP_CANONICAL_VARIABLES:
        # Check canonical name or name with _pm suffix
        h5_var_name = var if pass_upper == 'AM' else f"{var}_pm"
        if h5_var_name not in grp.variables and var in grp.variables:
            h5_var_name = var

        if h5_var_name in grp.variables:
            arr = grp.variables[h5_var_name][:]
            # Replace fill values (-9999.0) with NaN for floating point variables
            if np.issubdtype(arr.dtype, np.floating):
                fill_val = getattr(grp.variables[h5_var_name], '_FillValue', -9999.0)
                arr = np.where(arr == fill_val, np.nan, arr)
            data[var] = arr
        else:
            raise KeyError(f"Variable {h5_var_name} (canonical {var}) missing in {group_name}")

    return data


def extract_smap_raw(filepath: str, pass_name: str = 'AM') -> Dict[str, np.ndarray]:
    """Convenience wrapper to open and extract variables for an overpass."""
    ds = open_smap_h5(filepath)
    try:
        data = read_smap_pass_data(ds, pass_name)
    finally:
        ds.close()
    return data
