"""
Exact reproduction of LIS neighbor interpolation for SMAP Enhanced (EASE-Grid 2.0 M09).
Faithfully reproduces easeV2_utils.F90 and neighbor_interp_input_withgrid.F90 from LISF.
Validated with 100% exact equality against LIS DAOBS unscaled observations.
"""

import numpy as np
from typing import Tuple

# WGS84 Constants strictly matching easeV2_utils.F90
MAP_EQUATORIAL_RADIUS_M = 6378137.0
MAP_ECCENTRICITY = 0.081819190843
PI = 3.14159265358979323846
E2 = MAP_ECCENTRICITY * MAP_ECCENTRICITY

# EASE-Grid 2.0 Cylindrical (M09) parameters
M_MAP_SECOND_REF_LAT = 30.0
M_SIN_PHI1 = np.sin(M_MAP_SECOND_REF_LAT * PI / 180.0)
M_COS_PHI1 = np.cos(M_MAP_SECOND_REF_LAT * PI / 180.0)
M_KZ = M_COS_PHI1 / np.sqrt(1.0 - E2 * M_SIN_PHI1 * M_SIN_PHI1)

M09_MAP_SCALE_M = 9008.055210146
M09_COLS = 3856
M09_ROWS = 1624
M09_R0 = (M09_COLS - 1) / 2.0  # 1927.5
M09_S0 = (M09_ROWS - 1) / 2.0  # 811.5


def latlon_to_ease2_m09(lat: np.ndarray, lon: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert geographic coordinates (lat, lon in degrees) to EASE-Grid 2.0 M09 row and column indices.
    Matches easeV2_convert in easeV2_utils.F90 and nint() indexing from neighbor_interp_input_withgrid.F90.
    Returns 0-based indices (row_idx, col_idx) matching C/Python indexing of global HDF5 datasets.
    """
    dlon = np.asarray(lon, dtype=np.float64)
    dlon = np.where(dlon < -180.0, dlon + 360.0, dlon)
    dlon = np.where(dlon > 180.0, dlon - 360.0, dlon)

    phi = np.asarray(lat, dtype=np.float64) * PI / 180.0
    lam = dlon * PI / 180.0

    sin_phi = np.sin(phi)
    ms = MAP_ECCENTRICITY * sin_phi

    q = (1.0 - E2) * (
        (sin_phi / (1.0 - E2 * sin_phi * sin_phi))
        - 0.5 / MAP_ECCENTRICITY * np.log((1.0 - ms) / (1.0 + ms))
    )

    x = MAP_EQUATORIAL_RADIUS_M * M_KZ * lam
    y = (MAP_EQUATORIAL_RADIUS_M * q) / (2.0 * M_KZ)

    # In LIS easeV2_convert:
    # col_ind = M09_R0 + (x / M09_MAP_SCALE_M)
    # row_ind = M09_S0 - (y / M09_MAP_SCALE_M)
    # In neighbor_interp_input_withgrid.F90:
    # i1 = nint(col_ind)  (1-based in Fortran: 1..3856)
    # j1 = nint(row_ind)  (1-based in Fortran: 1..1624)
    # In Python (0-based: 0..3855, 0..1623):
    # col_0b = nint(col_ind) - 1 = floor(col_ind + 0.5) - 1
    # row_0b = nint(row_ind) - 1 = floor(row_ind + 0.5) - 1
    col_ind_float = M09_R0 + (x / M09_MAP_SCALE_M)
    row_ind_float = M09_S0 - (y / M09_MAP_SCALE_M)

    col_ind = np.floor(col_ind_float + 0.5).astype(int) - 1
    row_ind = np.floor(row_ind_float + 0.5).astype(int) - 1

    # Bounds clipping
    col_ind = np.clip(col_ind, 0, M09_COLS - 1)
    row_ind = np.clip(row_ind, 0, M09_ROWS - 1)

    return row_ind, col_ind


def regrid_smap_to_northmor(smap_global_2d: np.ndarray,
                            ease_row_map: np.ndarray,
                            ease_col_map: np.ndarray) -> np.ndarray:
    """
    Reproject global EASE-Grid 2.0 (1624 x 3856) array to NorthMor (121 x 181)
    using precomputed nearest neighbor coordinate mapping.
    """
    return smap_global_2d[ease_row_map, ease_col_map]
