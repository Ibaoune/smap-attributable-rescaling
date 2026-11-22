"""
SAR Quality Control (QC) and Masking Logic.
Implements SMAP recommended retrieval flag test and reproduces
all Noah-MP / LIS environmental and physical check filters (noahmp401_qc_soilmobs.F90).
"""

import numpy as np


def is_smap_recommended(retrieval_qual_flag: np.ndarray) -> np.ndarray:
    """
    Check if SMAP retrieval quality is recommended (bit 0 == 0).
    Returns boolean array where True indicates recommended quality.
    """
    return (retrieval_qual_flag.astype(np.uint16) & 1) == 0


def check_lis_vwc_exceeded(vwc: np.ndarray, threshold: float = 5.0) -> np.ndarray:
    """
    LIS SMAP plugin check (read_NASASMAPsm.F90:1155):
    Reject if VWC > 5.0 kg/m2.
    """
    return vwc > threshold


def check_lis_rain(rain_rate: np.ndarray, threshold: float = 3e-6) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if incident rain rate > 3e-6 kg/m2/s.
    """
    return rain_rate > threshold


def check_lis_frozen(smc: np.ndarray, sh2o: np.ndarray,
                     stc: np.ndarray, tkfrz: float = 273.15) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if soil is frozen (|smc - sh2o| > 1e-4 or soil temp <= 273.15 K).
    """
    ice_present = np.abs(smc - sh2o) > 1e-4
    temp_frozen = stc <= tkfrz
    return ice_present | temp_frozen


def check_lis_forest(vegt: np.ndarray) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if dominant land cover is forest (IGBP classes 1 to 4).
    """
    return (vegt >= 1) & (vegt <= 4)


def check_lis_snow(swe_m: np.ndarray, sca: np.ndarray) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if snow water equivalent > 1 mm (0.001 m) or snow cover area > 1e-4.
    """
    return (swe_m > 0.001) | (sca > 1e-4)


def check_lis_high_vegetation(shdfac: np.ndarray, threshold: float = 0.9) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if green vegetation fraction shdfac > 0.9.
    """
    return shdfac > threshold


def check_lis_tail_saturation(smcmax: np.ndarray, smobs: np.ndarray, margin: float = 0.02) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if observation is within margin of saturation (smcmax - smobs < 0.02).
    """
    return (smcmax - smobs) < margin


def check_lis_tail_wilting(smcwlt: np.ndarray, smobs: np.ndarray, margin: float = 0.02) -> np.ndarray:
    """
    LIS Noah-MP check: Reject if observation is within margin of wilting point (smobs - smcwlt < 0.02).
    """
    return (smobs - smcwlt) < margin
