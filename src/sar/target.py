"""
SAR Target Computation and Validation Module.
Defines d = SMAP_raw - SSM_OL (signed bias).
"""

import numpy as np


def compute_signed_target(smap_sm: np.ndarray, ssm_ol_daily: np.ndarray) -> np.ndarray:
    """
    Compute signed discrepancy:
    d = smap_sm - ssm_ol_daily
    Returns NaN where either observation or OL is missing.
    """
    diff = smap_sm - ssm_ol_daily
    # Mask invalid physical ranges if needed
    invalid = (smap_sm < 0.0) | (smap_sm > 1.0) | (ssm_ol_daily < 0.0) | (ssm_ol_daily > 1.0)
    diff = np.where(invalid, np.nan, diff)
    return diff

# Alias for compatibility
compute_target_discrepancy = compute_signed_target
