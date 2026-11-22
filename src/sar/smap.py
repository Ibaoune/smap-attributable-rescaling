"""
SMAP SPL3SMP_E Raw HDF5 Data Reader.
Handles extraction of AM and PM retrieval fields, quality flags, and geophysical variables.
Standardizes PM variable names (which carry the '_pm' suffix in HDF5) to canonical variable names.
Ensures masked arrays and fill values are properly converted to NaN with no residual -9999 or out-of-range values.
Supports optional bounding box slicing for orders-of-magnitude faster I/O.
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


def read_smap_pass_data(
    ds: netCDF4.Dataset, 
    pass_name: str = 'AM',
    bbox: Optional[Tuple[int, int, int, int]] = None
) -> Dict[str, np.ndarray]:
    """
    Extract all required SMAP variables for a specific overpass (AM or PM).
    Standardizes variable names (removes '_pm' suffix for PM pass).
    Converts masked arrays and fill values to NaN using np.ma.filled(arr.astype(float), np.nan).
    Enforces valid_min and valid_max if defined in the variable attributes.
    If bbox=(r_min, r_max, c_min, c_max) is provided, slices directly during HDF5 read.
    Returns a dictionary of 2D numpy arrays.
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
            var_obj = grp.variables[h5_var_name]
            if bbox is not None:
                r0, r1, c0, c1 = bbox
                raw_arr = var_obj[r0:r1+1, c0:c1+1]
            else:
                raw_arr = var_obj[:]
            
            # For floating point variables, convert masked array to float with NaN fill
            if np.issubdtype(var_obj.dtype, np.floating) or var in [
                'soil_moisture', 'vegetation_opacity', 'vegetation_water_content',
                'surface_temperature', 'roughness_coefficient', 'static_water_body_fraction'
            ]:
                if isinstance(raw_arr, np.ma.MaskedArray):
                    arr = np.ma.filled(raw_arr.astype(float), np.nan)
                else:
                    arr = raw_arr.astype(float)
                
                # Double-check fill value (-9999.0 or attribute)
                fill_val = getattr(var_obj, '_FillValue', -9999.0)
                arr[arr == fill_val] = np.nan
                arr[arr == -9999.0] = np.nan
                
                # Enforce valid range if present
                vmin = getattr(var_obj, 'valid_min', None)
                vmax = getattr(var_obj, 'valid_max', None)
                if vmin is not None:
                    arr[arr < vmin] = np.nan
                if vmax is not None:
                    arr[arr > vmax] = np.nan
            else:
                # Integer flags (retrieval_qual_flag, surface_flag)
                if isinstance(raw_arr, np.ma.MaskedArray):
                    arr = np.ma.filled(raw_arr, getattr(var_obj, '_FillValue', 65534))
                else:
                    arr = raw_arr
            data[var] = arr
        else:
            raise KeyError(f"Variable {h5_var_name} (canonical {var}) missing in {group_name}")

    return data


def extract_smap_raw(filepath: str, pass_name: str = 'AM', bbox: Optional[Tuple[int, int, int, int]] = None) -> Dict[str, np.ndarray]:
    """Convenience wrapper to open and extract variables for an overpass."""
    ds = open_smap_h5(filepath)
    try:
        data = read_smap_pass_data(ds, pass_name, bbox=bbox)
    finally:
        ds.close()
    return data
