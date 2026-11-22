"""
SAR Predictor Extraction and Feature Engineering Module.
Derives static features (continuous slope) and process features (rolling precip, dry spells,
relative moisture to SOILPARM.TBL limits, seasonal harmonics).
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple


def compute_topographic_slope(elevation: np.ndarray, res_deg: float = 0.05,
                              lat_ref: float = 33.0) -> np.ndarray:
    """
    Compute continuous topographic slope in degrees from 2D elevation grid.
    Uses central finite differences with metric scaling based on latitude.
    """
    dy_m = res_deg * 111139.0  # meters per degree latitude (~5557 m)
    dx_m = res_deg * 111139.0 * np.cos(np.radians(lat_ref))  # meters per degree longitude

    grad_y, grad_x = np.gradient(elevation, dy_m, dx_m)
    slope_rad = np.arctan(np.sqrt(grad_x**2 + grad_y**2))
    slope_deg = np.degrees(slope_rad)
    return slope_deg.astype(np.float32)


def compute_rolling_precip(precip_series: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """
    Compute antecedent precipitation rolling sums for 1d, 3d, and 7d.
    Note: strictly antecedent (lagged by 1 day) to prevent target leakage.
    """
    p_lag1 = precip_series.shift(1)
    p_3d = precip_series.shift(1).rolling(window=3, min_periods=1).sum()
    p_7d = precip_series.shift(1).rolling(window=7, min_periods=1).sum()
    return p_lag1, p_3d, p_7d


def compute_consecutive_dry_days(precip_series: pd.Series, threshold_mm: float = 1.0) -> pd.Series:
    """
    Calculate the number of consecutive days since the last rainfall event > threshold_mm.
    Strictly antecedent (shifts by 1 day) to avoid leakage.
    """
    is_rain = (precip_series.shift(1) > threshold_mm).fillna(False)
    dry_days = np.zeros(len(precip_series), dtype=np.int16)
    count = 0
    for i, rain in enumerate(is_rain):
        if rain:
            count = 0
        else:
            count += 1
        dry_days[i] = count
    return pd.Series(dry_days, index=precip_series.index)


def compute_relative_soil_moisture(ssm: np.ndarray, smcwlt: np.ndarray,
                                   smcref: np.ndarray) -> np.ndarray:
    """
    Calculate the relative position of SSM between wilting point and field capacity:
    rel_ssm = (ssm - smcwlt) / (smcref - smcwlt)
    Values are clipped to [0, 1.5] for robust modeling.
    """
    denom = smcref - smcwlt
    denom = np.where(denom <= 1e-4, np.nan, denom)
    rel = (ssm - smcwlt) / denom
    rel = np.clip(rel, 0.0, 1.5)
    return rel.astype(np.float32)


def compute_seasonal_harmonics(doy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute annual seasonal harmonics (sin and cos of DOY).
    """
    radians = 2.0 * np.pi * doy / 365.25
    return np.sin(radians).astype(np.float32), np.cos(radians).astype(np.float32)
