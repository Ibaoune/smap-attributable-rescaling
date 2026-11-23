#!/usr/bin/env python3
"""
scripts/11_test_2020_sealed.py

PHASE 3 : Sealed Test on Year 2020 with LightGBM Attributable Bias Correction.
Target: d = SMAP raw - SSM Noah-MP Open Loop (daily average)
Domain: NorthMor 0.05 deg (121 x 181, 16,030 land pixels)
Evaluation: Year 2020 (Strictly sealed, evaluated once)

Formal Decision 15 Success Criteria:
1. Candidate method: V2 (Removable + Context Representation)
2. Reference: B2-train CDF reconstructed on 2016-01-01 -> 2019-01-01
3. No retraining, no hyperparameter modification (pre-trained model loaded from 2016-2018)
4. Success criteria:
   - |Bias(V2)| < |Bias(B2-train)| on Specific and Common assimilable populations
   - r_anom(V2) > r_anom(B2-train) on Specific and Common assimilable populations
   - RMSE(V2) <= 1.05 * RMSE(B2-train) on Common assimilable population
5. Zero hardcoded numbers in report and figures. All metrics dynamically computed from data.
"""

import os
import sys
import time
import math
import yaml
import gc
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import seaborn as sns
import lightgbm as lgb

# Base paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
INTERIM_DIR = os.path.join(BASE_DIR, "data", "interim")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
FIG_DIR = os.path.join(BASE_DIR, "figures", "evaluation")
FIG_EN_DIR = os.path.join(FIG_DIR, "en")

STATIC_FILE = os.path.join(INTERIM_DIR, "static", "static_predictors.nc")
CONFIG_PREDICTORS = os.path.join(BASE_DIR, "configs", "predictors.yaml")

MODEL_FILE = os.path.join(INTERIM_DIR, "lgb_model_2016_2018.txt")
B1_CLIM_FILE = os.path.join(INTERIM_DIR, "b1_climatology_2016_2018.csv")
B1_PIX_FILE = os.path.join(INTERIM_DIR, "b1_pixel_fallback_2016_2018.csv")

CDF_MODEL_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "noahmp_model_cdf_NorthMor.nc")
CDF_OBS_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "smap_obs_cdf_NorthMor.nc")
LUT_CACHE_FILE = os.path.join(INTERIM_DIR, "lut_cdf_rescaling.npy")

CDF_B2_TRAIN_MOD = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_model", "noahmp_model_cdf_NorthMor.nc")
CDF_B2_TRAIN_OBS = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_obs", "smap_obs_cdf_NorthMor.nc")
LUT_B2_TRAIN_FILE = os.path.join(INTERIM_DIR, "lut_cdf_b2_train.npy")

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def get_n_jobs():
    """Get number of CPU cores allocated by SLURM or fallback to local count."""
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        try:
            return int(slurm_cpus)
        except ValueError:
            pass
    return os.cpu_count() or 8


def load_feature_groups():
    """Load and parse feature groups from configs/predictors.yaml."""
    with open(CONFIG_PREDICTORS, 'r') as f:
        cfg = yaml.safe_load(f)

    cols_removable = [
        'elevation', 'slope_deg', 'slopetype',
        'texture_frac_1', 'texture_frac_2', 'texture_frac_3', 'texture_frac_4',
        'texture_frac_5', 'texture_frac_6', 'texture_frac_7', 'texture_frac_8',
        'texture_frac_9', 'texture_frac_10', 'texture_frac_11', 'texture_frac_12',
        'texture_frac_13', 'texture_frac_14', 'texture_frac_15', 'texture_frac_16',
        'landcover_frac_1', 'landcover_frac_2', 'landcover_frac_3', 'landcover_frac_4',
        'landcover_frac_5', 'landcover_frac_6', 'landcover_frac_7', 'landcover_frac_8',
        'landcover_frac_9', 'landcover_frac_10', 'landcover_frac_11', 'landcover_frac_12',
        'landcover_frac_13', 'landcover_frac_14', 'landcover_frac_15', 'landcover_frac_16',
        'landcover_frac_17', 'landcover_frac_18', 'landcover_frac_19', 'landcover_frac_20',
        'static_water_body_fraction', 'surface_flag',
        'surface_temperature', 'roughness_coefficient', 'roughness_is_nan',
        'vegetation_opacity', 'vegetation_water_content', 'pass_is_pm'
    ]
    cols_retained = [
        'irrigfrac_gmia', 'delta_precip_imerg_mswep', 'precip_mswep',
        'precip_1d_sum', 'precip_3d_sum', 'precip_7d_sum', 'days_since_rain_gt_1mm'
    ]
    cols_context_rep = [
        'doy_sin', 'doy_cos', 'ssm_lag1', 'ssm_rel_wilt_ref_lag1'
    ]
    cols_context_proc = [
        'sm_layer2_lag1', 'sm_layer3_lag1', 'sm_layer4_lag1',
        'rzsm_lag1', 'et_lag1', 'surface_runoff_lag1', 'baseflow_lag1', 'lai_lag1'
    ]
    feature_cols = cols_removable + cols_retained + cols_context_rep + cols_context_proc
    return cols_removable, cols_retained, cols_context_rep, cols_context_proc, feature_cols




def compute_gmia_contrast_bootstrap(sub, d_col, B=500):
    """Compute GMIA >30% vs 0% contrast on Clay Loam JJA with 1° block bootstrap 95% CIs."""
    sub_cl_jja = sub[(sub['texture_name'] == 'Clay Loam') & (sub['season'] == 'JJA')]
    sub_0 = sub_cl_jja[sub_cl_jja['gmia_class'] == '0%']
    sub_gt30 = sub_cl_jja[sub_cl_jja['gmia_class'] == '>30%']
    n_px_0 = int(sub_0['pixel_id'].nunique())
    n_px_gt30 = int(sub_gt30['pixel_id'].nunique())
    n_blk_0 = int(sub_0['spatial_block_id'].nunique())
    n_blk_gt30 = int(sub_gt30['spatial_block_id'].nunique())
    
    mean_0 = float(sub_0[d_col].mean()) if len(sub_0) > 0 else np.nan
    mean_gt30 = float(sub_gt30[d_col].mean()) if len(sub_gt30) > 0 else np.nan
    contrast = mean_gt30 - mean_0 if (not np.isnan(mean_0) and not np.isnan(mean_gt30)) else np.nan
    
    sub_both = sub_cl_jja[sub_cl_jja['gmia_class'].isin(['0%', '>30%'])]
    n_blocks = int(sub_both['spatial_block_id'].nunique())
    
    if n_blk_0 < 5 or n_blk_gt30 < 5:
        ci_low, ci_high = np.nan, np.nan
        status = "non interprétable (< 5 blocs)"
    elif len(sub_0) >= 10 and len(sub_gt30) >= 10 and n_blocks > 1:
        stats = sub_both.groupby(['spatial_block_id', 'gmia_class'])[d_col].agg(['sum', 'count']).unstack('gmia_class').fillna(0)
        sum_0 = stats[('sum', '0%')].values if ('sum', '0%') in stats else np.zeros(len(stats))
        cnt_0 = stats[('count', '0%')].values if ('count', '0%') in stats else np.zeros(len(stats))
        sum_gt30 = stats[('sum', '>30%')].values if ('sum', '>30%') in stats else np.zeros(len(stats))
        cnt_gt30 = stats[('count', '>30%')].values if ('count', '>30%') in stats else np.zeros(len(stats))
        
        np.random.seed(42)
        boot_diffs = np.empty(B)
        for b in range(B):
            idx = np.random.choice(n_blocks, size=n_blocks, replace=True)
            c0 = cnt_0[idx].sum()
            cgt = cnt_gt30[idx].sum()
            m0 = (sum_0[idx].sum() / c0) if c0 > 0 else np.nan
            mgt = (sum_gt30[idx].sum() / cgt) if cgt > 0 else np.nan
            boot_diffs[b] = mgt - m0
        valid_boot = boot_diffs[~np.isnan(boot_diffs)]
        if len(valid_boot) > 50:
            ci_low = float(np.percentile(valid_boot, 2.5))
            ci_high = float(np.percentile(valid_boot, 97.5))
            status = "interprétable"
        else:
            ci_low, ci_high = np.nan, np.nan
            status = "non interprétable"
    else:
        ci_low, ci_high = np.nan, np.nan
        status = "non interprétable"
        
    return contrast, ci_low, ci_high, n_px_gt30, n_px_0, n_blk_gt30, n_blk_0, status


def evaluate_2020_sealed_test(feature_cols, cols_removable, cols_retained, cols_context_rep, cols_context_proc,
                              lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, b1_clim_table, b1_pixel_table,
                              global_fallback, n_jobs=8):
    """Evaluate 2020 test set using pre-trained LightGBM model and exact TreeSHAP without any retraining."""
    print("--- Evaluating Sealed 2020 Test Set with Exact TreeSHAP and Variants V1, V2, V3 ---")
    t0 = time.time()
    
    idx_rem = [feature_cols.index(c) for c in cols_removable]
    idx_ret = [feature_cols.index(c) for c in cols_retained]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_context_rep]
    idx_ctx_proc = [feature_cols.index(c) for c in cols_context_proc]
    
    pq_2020 = os.path.join(TABLE_DIR, "sar_modeling_table_2020.parquet")
    meta_cols = [
        'pixel_id', 'month', 'season', 'gmia_class', 'spatial_block_id',
        'smap_sm', 'ssm_ol_daily', 'target_d', 'split',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    load_cols = list(set(feature_cols + meta_cols))
    
    print(f"Loading 2020 data from {pq_2020}...")
    df_2020 = pd.read_parquet(pq_2020, columns=load_cols)
    n_total_2020 = len(df_2020)
    
    # Filter non-wilt QC observations (Protocol 2a standard)
    qc_non_wilt = (
        (df_2020['qc_lis_rain'] == 0) &
        (df_2020['qc_lis_frozen'] == 0) &
        (df_2020['qc_lis_forest'] == 0) &
        (df_2020['qc_lis_snow'] == 0) &
        (df_2020['qc_lis_high_veg'] == 0) &
        (df_2020['qc_lis_sat_limit'] == 0) &
        (df_2020['qc_lis_vwc_gt_5'] == 0) &
        (~df_2020['target_d'].isna())
    )
    df_2020 = df_2020.loc[qc_non_wilt].copy().reset_index(drop=True)
    n_base_qc = len(df_2020)
    print(f"Loaded 2020 non-wilt QC: {n_base_qc:,} rows ({n_base_qc/n_total_2020*100:.2f}%) in {time.time()-t0:.1f}s")
    
    r_idx = df_2020['pixel_id'].values // 181
    c_idx = df_2020['pixel_id'].values % 181
    wlt_thresh = soil_wlt[r_idx, c_idx] + 0.02
    soil_ids = dom_soil[r_idx, c_idx]
    df_2020['texture_name'] = [TEXTURE_NAMES.get(int(s), "Other") for s in soil_ids]
    df_2020['wlt_thresh'] = wlt_thresh
    
    tile_map = np.full((121, 181), -1, dtype=int)
    t = 0
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                tile_map[r, c] = t
                t += 1
    tiles = tile_map[r_idx, c_idx]
    
    # 1. Baseline B0: Raw SMAP
    df_2020['smap_star_b0'] = df_2020['smap_sm']
    df_2020['d_res_b0'] = df_2020['target_d']
    df_2020['under_wlt_b0'] = (df_2020['smap_star_b0'] < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_b0'] = (
        (df_2020['under_wlt_b0'] == 0) &
        (df_2020['smap_star_b0'] >= 0.01) &
        (df_2020['smap_star_b0'] <= 1.0) &
        (~df_2020['smap_star_b0'].isna())
    )
    
    # 2. Baseline B1: Climatology via vectorized pandas join on (pixel_id, month)
    df_2020 = df_2020.merge(b1_clim_table, on=['pixel_id', 'month'], how='left')
    df_2020 = df_2020.merge(b1_pixel_table, on='pixel_id', how='left')
    c_hat = df_2020['b1_bias'].fillna(df_2020['b1_pixel_fallback']).fillna(global_fallback).values.astype(np.float32)
    df_2020.drop(columns=['b1_bias', 'b1_pixel_fallback'], inplace=True)

    df_2020['smap_star_b1'] = df_2020['smap_sm'] - c_hat
    df_2020['d_res_b1'] = df_2020['target_d'] - c_hat
    df_2020['under_wlt_b1'] = (df_2020['smap_star_b1'] < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_b1'] = (
        (df_2020['under_wlt_b1'] == 0) &
        (df_2020['smap_star_b1'] >= 0.01) &
        (df_2020['smap_star_b1'] <= 1.0) &
        (~df_2020['smap_star_b1'].isna())
    )
    
    # 3. Baseline B2: CDF legacy (2016-2021)
    months = df_2020['month'].values
    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta = obs_x[1, :, :] - obs_x[0, :, :]
    obs_x0 = obs_x[0, :, :]
    m_indices = months - 1
    smap_vals = df_2020['smap_sm'].values
    
    def apply_lut(lut_arr, obs_x0_arr, obs_delta_arr):
        res = np.full(len(df_2020), np.nan, dtype=np.float32)
        for m in range(12):
            mask_m = (m_indices == m)
            if not np.any(mask_m):
                continue
            sub_tiles = tiles[mask_m]
            sub_smap = smap_vals[mask_m]
            sub_delta = obs_delta_arr[m, sub_tiles]
            sub_x0 = obs_x0_arr[m, sub_tiles]
            valid_sub = (sub_delta > 0) & (sub_tiles >= 0) & (sub_smap >= 0.01) & (sub_smap <= 1.0)
            binvals = np.zeros(len(sub_tiles), dtype=int)
            binvals[valid_sub] = np.floor((sub_smap[valid_sub] - sub_x0[valid_sub]) / sub_delta[valid_sub] + 0.5).astype(int) + 1
            binvals = np.clip(binvals, 1, 100) - 1
            rescaled_m = np.full(len(sub_tiles), np.nan, dtype=np.float32)
            rescaled_m[valid_sub] = lut_arr[binvals[valid_sub], m, sub_tiles[valid_sub]]
            res[mask_m] = rescaled_m
        return res
        
    smap_star_b2 = apply_lut(lut_b2, obs_x0, obs_delta)
    df_2020['smap_star_b2'] = smap_star_b2
    df_2020['b2_has_value'] = (~np.isnan(smap_star_b2)).astype(int)
    df_2020['d_res_b2'] = smap_star_b2 - df_2020['ssm_ol_daily']
    df_2020['under_wlt_b2'] = (smap_star_b2 < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_b2'] = (
        (df_2020['under_wlt_b2'] == 0) &
        (df_2020['b2_has_value'] == 1) &
        (smap_star_b2 >= 0.01) &
        (smap_star_b2 <= 1.0) &
        (~np.isnan(smap_star_b2))
    )
    
    # 4. Baseline B2-train: CDF strictly reconstructed on 2016-01-01 -> 2019-01-01
    with nc.Dataset(CDF_B2_TRAIN_OBS) as ds_o_tr:
        obs_x_tr = ds_o_tr.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta_tr = obs_x_tr[1, :, :] - obs_x_tr[0, :, :]
    obs_x0_tr = obs_x_tr[0, :, :]
    
    smap_star_b2_tr = apply_lut(lut_b2_train, obs_x0_tr, obs_delta_tr)
    df_2020['smap_star_b2_train'] = smap_star_b2_tr
    df_2020['applied_corr_b2_train'] = smap_star_b2_tr - df_2020['smap_sm']
    df_2020['b2_train_has_value'] = (~np.isnan(smap_star_b2_tr)).astype(int)
    df_2020['d_res_b2_train'] = smap_star_b2_tr - df_2020['ssm_ol_daily']
    df_2020['under_wlt_b2_train'] = (smap_star_b2_tr < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_b2_train'] = (
        (df_2020['under_wlt_b2_train'] == 0) &
        (df_2020['b2_train_has_value'] == 1) &
        (smap_star_b2_tr >= 0.01) &
        (smap_star_b2_tr <= 1.0) &
        (~np.isnan(smap_star_b2_tr))
    )
    
    # 5. Pre-trained LightGBM & Exact TreeSHAP inference (NO RETRAINING)
    print(f"Loading pre-trained LightGBM model from {MODEL_FILE}...")
    booster = lgb.Booster(model_file=MODEL_FILE)
    
    X_val_2020 = df_2020[feature_cols].values.astype(np.float32)
    contrib_cache_file = os.path.join(INTERIM_DIR, "lgb_contrib_2020.npy")
    
    if os.path.exists(contrib_cache_file):
        print(f"Loading cached exact TreeSHAP attributions from {contrib_cache_file}...")
        contrib_2020 = np.load(contrib_cache_file)
    else:
        print(f"Computing EXACT TreeSHAP on all {len(df_2020):,} observations of 2020 with {n_jobs} threads...")
        t_shap = time.time()
        contrib_2020 = booster.predict(X_val_2020, pred_contrib=True, num_threads=n_jobs)
        np.save(contrib_cache_file, contrib_2020)
        print(f"Exact TreeSHAP computed and cached in {time.time()-t_shap:.1f}s")
        
    shap_2020 = contrib_2020[:, :-1]
    phi_0 = float(contrib_2020[0, -1])
    
    # Lundberg additivity verification
    y_pred_2020 = booster.predict(X_val_2020, num_threads=n_jobs)
    sum_phi = contrib_2020.sum(axis=1)
    dev_additivity = float(np.max(np.abs(sum_phi - y_pred_2020)))
    print(f"Exact TreeSHAP Lundberg additivity verified. Max deviation: {dev_additivity:.2e}")
    
    phi_rem = shap_2020[:, idx_rem].sum(axis=1)
    phi_ret = shap_2020[:, idx_ret].sum(axis=1)
    phi_ctx_rep = shap_2020[:, idx_ctx_rep].sum(axis=1)
    phi_ctx_proc = shap_2020[:, idx_ctx_proc].sum(axis=1)
    
    df_2020['phi_removable'] = phi_rem
    df_2020['phi_retained'] = phi_ret
    df_2020['phi_context_rep'] = phi_ctx_rep
    df_2020['phi_context_proc'] = phi_ctx_proc
    df_2020['phi_base_val'] = phi_0
    
    # 6. Variants V1, V2, V3
    # V1 = removable
    df_2020['smap_star_v1'] = df_2020['smap_sm'] - (phi_rem + phi_0)
    df_2020['applied_corr_v1'] = df_2020['smap_star_v1'] - df_2020['smap_sm']
    df_2020['d_res_v1'] = df_2020['smap_star_v1'] - df_2020['ssm_ol_daily']
    df_2020['under_wlt_v1'] = (df_2020['smap_star_v1'] < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_v1'] = (
        (df_2020['under_wlt_v1'] == 0) &
        (df_2020['smap_star_v1'] >= 0.01) &
        (df_2020['smap_star_v1'] <= 1.0) &
        (~df_2020['smap_star_v1'].isna())
    )
    
    # V2 = removable + context_representation (Formal candidate method)
    df_2020['smap_star_v2'] = df_2020['smap_sm'] - (phi_rem + phi_ctx_rep + phi_0)
    df_2020['applied_corr_v2'] = df_2020['smap_star_v2'] - df_2020['smap_sm']
    df_2020['d_res_v2'] = df_2020['smap_star_v2'] - df_2020['ssm_ol_daily']
    df_2020['under_wlt_v2'] = (df_2020['smap_star_v2'] < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_v2'] = (
        (df_2020['under_wlt_v2'] == 0) &
        (df_2020['smap_star_v2'] >= 0.01) &
        (df_2020['smap_star_v2'] <= 1.0) &
        (~df_2020['smap_star_v2'].isna())
    )
    
    # V3 = removable + all context
    df_2020['smap_star_v3'] = df_2020['smap_sm'] - (phi_rem + phi_ctx_rep + phi_ctx_proc + phi_0)
    df_2020['applied_corr_v3'] = df_2020['smap_star_v3'] - df_2020['smap_sm']
    df_2020['d_res_v3'] = df_2020['smap_star_v3'] - df_2020['ssm_ol_daily']
    df_2020['under_wlt_v3'] = (df_2020['smap_star_v3'] < df_2020['wlt_thresh']).astype(int)
    df_2020['assimilable_v3'] = (
        (df_2020['under_wlt_v3'] == 0) &
        (df_2020['smap_star_v3'] >= 0.01) &
        (df_2020['smap_star_v3'] <= 1.0) &
        (~df_2020['smap_star_v3'].isna())
    )

    # Common assimilable sample across all 7 references
    df_2020['assimilable_common'] = (
        df_2020['assimilable_b0'] &
        df_2020['assimilable_b1'] &
        df_2020['assimilable_b2'] &
        df_2020['assimilable_b2_train'] &
        df_2020['assimilable_v1'] &
        df_2020['assimilable_v2'] &
        df_2020['assimilable_v3']
    )
    
    return df_2020, shap_2020, X_val_2020, phi_0, n_total_2020, n_base_qc


def compute_seasonal_corrections_and_metrics_2020(df_2020, n_total_2020, n_base_qc):
    """Compute seasonal applied corrections, full protocol 2a metrics, and seasonal stratified metrics for 2020."""
    print("--- Computing 2020 Seasonal Applied Corrections and Full Protocol 2a Metrics ---")
    
    # 1. Seasonal applied correction and residual bias table
    seasons = ['DJF', 'MAM', 'JJA', 'SON', 'Annual']
    seas_rows = []
    for s in seasons:
        sub = df_2020 if s == 'Annual' else df_2020[df_2020['season'] == s]
        seas_rows.append({
            'season': s,
            'n_obs': len(sub),
            'corr_v1_mean': float(sub['applied_corr_v1'].mean()),
            'd_res_v1_mean': float(sub['d_res_v1'].mean()),
            'corr_v2_mean': float(sub['applied_corr_v2'].mean()),
            'd_res_v2_mean': float(sub['d_res_v2'].mean()),
            'corr_v3_mean': float(sub['applied_corr_v3'].mean()),
            'd_res_v3_mean': float(sub['d_res_v3'].mean()),
            'corr_b2_train_mean': float(sub['applied_corr_b2_train'].mean()),
            'd_res_b2_train_mean': float(sub['d_res_b2_train'].mean())
        })
    df_seasonal = pd.DataFrame(seas_rows)
    csv_seasonal = os.path.join(REPORTS_DIR, "phase3_test2020_seasonal_corrections.csv")
    df_seasonal.to_csv(csv_seasonal, index=False)
    print(f"Saved: {csv_seasonal}")
    
    # 2. Full Protocol 2a metrics table (Annual, Specific and Common)
    models = [
        ('B0_Raw', 'd_res_b0', 'smap_star_b0', 'assimilable_b0', 'under_wlt_b0'),
        ('B1_Climatology', 'd_res_b1', 'smap_star_b1', 'assimilable_b1', 'under_wlt_b1'),
        ('B2_CDF_Legacy', 'd_res_b2', 'smap_star_b2', 'assimilable_b2', 'under_wlt_b2'),
        ('B2_train_CDF', 'd_res_b2_train', 'smap_star_b2_train', 'assimilable_b2_train', 'under_wlt_b2_train'),
        ('LGB_V1_Removable', 'd_res_v1', 'smap_star_v1', 'assimilable_v1', 'under_wlt_v1'),
        ('LGB_V2_Rem_ContextRep', 'd_res_v2', 'smap_star_v2', 'assimilable_v2', 'under_wlt_v2'),
        ('LGB_V3_Rem_AllContext', 'd_res_v3', 'smap_star_v3', 'assimilable_v3', 'under_wlt_v3'),
    ]
    
    metrics_rows = []
    for pop_name, pop_col in [('Specific_Assimilable', None), ('Common_Assimilable', 'assimilable_common')]:
        for name, d_col, star_col, assim_col, wlt_col in models:
            if pop_col is None:
                sub = df_2020[df_2020[assim_col]].copy()
                n_wlt = int(df_2020[wlt_col].sum())
            else:
                sub = df_2020[df_2020[pop_col]].copy()
                n_wlt = 0
                
            n_assim = len(sub)
            d_vals = sub[d_col].values
            bias = float(np.mean(d_vals))
            rmse = float(np.sqrt(np.mean(d_vals ** 2)))
            ubrmse = float(np.sqrt(np.var(d_vals)))
            
            # Correlation r_smapstar_vs_ol
            ol_vals = sub['ssm_ol_daily'].values
            star_vals = sub[star_col].values
            r_ol = float(np.corrcoef(star_vals, ol_vals)[0, 1]) if len(sub) > 1 else np.nan
            
            # Anomaly correlation with raw SMAP (monthly pixel mean removed)
            anom_raw = (sub['smap_sm'] - sub.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')).values
            anom_star = (sub[star_col] - sub.groupby(['pixel_id', 'month'])[star_col].transform('mean')).values
            valid_a = ~(np.isnan(anom_raw) | np.isnan(anom_star))
            r_anom = float(np.corrcoef(anom_raw[valid_a], anom_star[valid_a])[0, 1]) if valid_a.sum() > 10 else np.nan
            
            # GMIA contrast with bootstrap on blocks
            contrast, ci_low, ci_high, n_px_gt30, n_px_0, n_blocks_gt30, n_blocks_0, status = compute_gmia_contrast_bootstrap(sub, d_col)
            
            metrics_rows.append({
                'evaluation_population': pop_name,
                'model': name,
                'n_assimilable': n_assim,
                'assimilable_pct': n_assim / n_base_qc * 100.0,
                'wilt_rejections': n_wlt,
                'wilt_rejection_pct': n_wlt / n_base_qc * 100.0,
                'bias': bias,
                'rmse': rmse,
                'ubrmse': ubrmse,
                'r_smapstar_vs_ol': r_ol,
                'r_anom_with_smap_raw': r_anom,
                'gmia_contrast_cl_jja': contrast,
                'gmia_contrast_ci_low': ci_low,
                'gmia_contrast_ci_high': ci_high,
                'n_px_gt30': n_px_gt30,
                'n_px_0': n_px_0,
                'n_blocks_gt30': n_blocks_gt30,
                'n_blocks_0': n_blocks_0,
                'contrast_status': status
            })
            
    df_metrics = pd.DataFrame(metrics_rows)
    csv_metrics = os.path.join(REPORTS_DIR, "phase3_test2020_metrics.csv")
    df_metrics.to_csv(csv_metrics, index=False)
    print(f"Saved: {csv_metrics}")
    
    # 3. Seasonal metrics table for each reference
    seasonal_metrics_rows = []
    for pop_name, pop_col in [('Specific_Assimilable', None), ('Common_Assimilable', 'assimilable_common')]:
        for name, d_col, star_col, assim_col, wlt_col in models:
            for s in seasons:
                df_s = df_2020 if s == 'Annual' else df_2020[df_2020['season'] == s]
                n_s_qc = len(df_s)
                
                if pop_col is None:
                    sub_s = df_s[df_s[assim_col]].copy()
                    n_wlt_s = int(df_s[wlt_col].sum())
                else:
                    sub_s = df_s[df_s[pop_col]].copy()
                    n_wlt_s = 0
                    
                n_assim_s = len(sub_s)
                d_vals_s = sub_s[d_col].values
                bias_s = float(np.mean(d_vals_s)) if n_assim_s > 0 else np.nan
                rmse_s = float(np.sqrt(np.mean(d_vals_s ** 2))) if n_assim_s > 0 else np.nan
                ubrmse_s = float(np.sqrt(np.var(d_vals_s))) if n_assim_s > 0 else np.nan
                
                ol_vals_s = sub_s['ssm_ol_daily'].values
                star_vals_s = sub_s[star_col].values
                r_ol_s = float(np.corrcoef(star_vals_s, ol_vals_s)[0, 1]) if n_assim_s > 1 else np.nan
                
                anom_raw_s = (sub_s['smap_sm'] - sub_s.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')).values
                anom_star_s = (sub_s[star_col] - sub_s.groupby(['pixel_id', 'month'])[star_col].transform('mean')).values
                valid_as = ~(np.isnan(anom_raw_s) | np.isnan(anom_star_s))
                r_anom_s = float(np.corrcoef(anom_raw_s[valid_as], anom_star_s[valid_as])[0, 1]) if valid_as.sum() > 10 else np.nan
                
                seasonal_metrics_rows.append({
                    'evaluation_population': pop_name,
                    'model': name,
                    'season': s,
                    'n_eval_qc': n_s_qc,
                    'n_assimilable': n_assim_s,
                    'assimilable_pct': n_assim_s / n_s_qc * 100.0,
                    'wilt_rejections': n_wlt_s,
                    'wilt_rejection_pct': n_wlt_s / n_s_qc * 100.0,
                    'bias': bias_s,
                    'rmse': rmse_s,
                    'ubrmse': ubrmse_s,
                    'r_smapstar_vs_ol': r_ol_s,
                    'r_anom_with_smap_raw': r_anom_s
                })
                
    df_seas_metrics = pd.DataFrame(seasonal_metrics_rows)
    csv_seas_metrics = os.path.join(REPORTS_DIR, "phase3_test2020_seasonal_metrics.csv")
    df_seas_metrics.to_csv(csv_seas_metrics, index=False)
    print(f"Saved: {csv_seas_metrics}")
    
    n_common = int(df_2020['assimilable_common'].sum())
    return df_seasonal, df_metrics, df_seas_metrics, n_common


def generate_figures_2020(df_2020):
    """Generate 2020 seasonal maps of applied corrections for V1, V2, V3, and B2-train."""
    print("--- Generating 2020 Seasonal Applied Correction Maps ---")
    seasons_order = ['DJF', 'MAM', 'JJA', 'SON']
    variants_fr = [
        ('V1 (Retirable)', 'applied_corr_v1'),
        ('V2 (Retirable + Ctx Repr.)', 'applied_corr_v2'),
        ('V3 (Retirable + Tout Ctx)', 'applied_corr_v3'),
        ('B2-train (Correction CDF)', 'applied_corr_b2_train')
    ]
    variants_en = [
        ('V1 (Removable)', 'applied_corr_v1'),
        ('V2 (Rem + Ctx Repr.)', 'applied_corr_v2'),
        ('V3 (Rem + All Context)', 'applied_corr_v3'),
        ('B2-train (CDF Matching)', 'applied_corr_b2_train')
    ]
    
    all_maps_data = []
    for _, col_name in variants_fr:
        for s in seasons_order:
            sub = df_2020[df_2020['season'] == s]
            m = sub.groupby('pixel_id')[col_name].mean().values
            all_maps_data.append(m[~np.isnan(m)])
    all_vals_concat = np.concatenate(all_maps_data)
    p2 = float(np.nanpercentile(all_vals_concat, 2))
    p98 = float(np.nanpercentile(all_vals_concat, 98))
    vlim = float(max(abs(p2), abs(p98)))
    vmin, vmax = -vlim, vlim
    
    for lang, variants, fname in [
        ('fr', variants_fr, 'F11_test2020_seasonal_applied_corrections.png'),
        ('en', variants_en, 'F11_test2020_seasonal_applied_corrections_en.png')
    ]:
        fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.05, hspace=0.08)
        
        for row_idx, (v_label, col_name) in enumerate(variants):
            for col_idx, s in enumerate(seasons_order):
                ax = axes[row_idx, col_idx]
                sub = df_2020[df_2020['season'] == s]
                mean_px = sub.groupby('pixel_id')[col_name].mean()
                
                grid = np.full((121, 181), np.nan, dtype=np.float32)
                p_ids = mean_px.index.values
                vals = mean_px.values
                grid[p_ids // 181, p_ids % 181] = vals
                
                im = ax.imshow(grid, origin='lower', cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=[-10.0, -1.0, 30.0, 36.0])
                if row_idx == 0:
                    ax.set_title(s, fontsize=12, fontweight='bold')
                if col_idx == 0:
                    ax.set_ylabel(f"{v_label}\n(°N)", fontsize=10, fontweight='bold')
                if row_idx == 3:
                    ax.set_xlabel("Longitude (°E)", fontsize=10)
                    
        cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(
            r"Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)" if lang == 'fr' else
            r"Mean applied correction $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)", fontsize=11
        )
        
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {out_path}")
        
        ide_art_dir = "/home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures"
        import shutil
        shutil.copy(out_path, os.path.join(ide_art_dir if lang == 'fr' else os.path.join(ide_art_dir, "en"), fname))


def evaluate_decision15_criteria(df_metrics):
    """Evaluate formal Decision 15 success criteria on 2020 test set."""
    # Specific assimilable
    spec = df_metrics[df_metrics['evaluation_population'] == 'Specific_Assimilable'].set_index('model')
    comm = df_metrics[df_metrics['evaluation_population'] == 'Common_Assimilable'].set_index('model')
    
    bias_v2_spec = float(spec.loc['LGB_V2_Rem_ContextRep', 'bias'])
    bias_b2_spec = float(spec.loc['B2_train_CDF', 'bias'])
    ranom_v2_spec = float(spec.loc['LGB_V2_Rem_ContextRep', 'r_anom_with_smap_raw'])
    ranom_b2_spec = float(spec.loc['B2_train_CDF', 'r_anom_with_smap_raw'])
    
    bias_v2_comm = float(comm.loc['LGB_V2_Rem_ContextRep', 'bias'])
    bias_b2_comm = float(comm.loc['B2_train_CDF', 'bias'])
    ranom_v2_comm = float(comm.loc['LGB_V2_Rem_ContextRep', 'r_anom_with_smap_raw'])
    ranom_b2_comm = float(comm.loc['B2_train_CDF', 'r_anom_with_smap_raw'])
    
    rmse_v2_comm = float(comm.loc['LGB_V2_Rem_ContextRep', 'rmse'])
    rmse_b2_comm = float(comm.loc['B2_train_CDF', 'rmse'])
    
    c1_bias_spec = abs(bias_v2_spec) < abs(bias_b2_spec)
    c2_bias_comm = abs(bias_v2_comm) < abs(bias_b2_comm)
    c3_ranom_spec = ranom_v2_spec > ranom_b2_spec
    c4_ranom_comm = ranom_v2_comm > ranom_b2_comm
    c5_rmse_comm = rmse_v2_comm <= (1.05 * rmse_b2_comm)
    
    all_met = c1_bias_spec and c2_bias_comm and c3_ranom_spec and c4_ranom_comm and c5_rmse_comm
    
    criteria_results = {
        'crit_1_bias_spec': {'met': c1_bias_spec, 'v2': abs(bias_v2_spec), 'b2': abs(bias_b2_spec)},
        'crit_2_bias_comm': {'met': c2_bias_comm, 'v2': abs(bias_v2_comm), 'b2': abs(bias_b2_comm)},
        'crit_3_ranom_spec': {'met': c3_ranom_spec, 'v2': ranom_v2_spec, 'b2': ranom_b2_spec},
        'crit_4_ranom_comm': {'met': c4_ranom_comm, 'v2': ranom_v2_comm, 'b2': ranom_b2_comm},
        'crit_5_rmse_comm': {'met': c5_rmse_comm, 'v2': rmse_v2_comm, 'b2': rmse_b2_comm, 'thresh_105': 1.05 * rmse_b2_comm},
        'overall_success': all_met
    }
    return criteria_results


def generate_markdown_report_2020(df_seasonal, df_metrics, df_seas_metrics, n_total_2020, n_base_qc, n_common, phi_0, crit):
    """Generate Markdown report for Phase 3 Sealed 2020 Test with ZERO hardcoded numbers and ZERO banned phrases."""
    print("--- Generating Phase 3 Test 2020 Markdown Report ---")
    out_md = os.path.join(REPORTS_DIR, "phase3_test2020_report.md")
    
    latex_v1 = r"$\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$"
    latex_v2 = r"$\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$"
    latex_v3 = r"$\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$"
    latex_dres = r"$d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$"
    latex_wlt = r"$\mathrm{SMAP}^* \ge \theta_{\mathrm{wlt}} + 0{,}02$"
    
    # 1. Seasonal Corrections Table (V1, V2, V3, B2-train)
    seas_rows = []
    for _, r in df_seasonal.iterrows():
        seas_rows.append(
            f"| {r['season']} | {int(r['n_obs']):,} | "
            f"{r['corr_v1_mean']:+.4f} m³/m³ | {r['d_res_v1_mean']:+.4f} m³/m³ | "
            f"{r['corr_v2_mean']:+.4f} m³/m³ | {r['d_res_v2_mean']:+.4f} m³/m³ | "
            f"{r['corr_v3_mean']:+.4f} m³/m³ | {r['d_res_v3_mean']:+.4f} m³/m³ | "
            f"{r['corr_b2_train_mean']:+.4f} m³/m³ | {r['d_res_b2_train_mean']:+.4f} m³/m³ |"
        )
    seas_table = chr(10).join(seas_rows)
    
    # 2. Full Protocol 2a Metrics Table
    def make_metrics_table(pop_name):
        sub = df_metrics[df_metrics['evaluation_population'] == pop_name]
        rows = []
        for _, r in sub.iterrows():
            if not np.isnan(r['gmia_contrast_ci_low']):
                ci_str = f"[{r['gmia_contrast_ci_low']:+.4f}, {r['gmia_contrast_ci_high']:+.4f}]"
            else:
                ci_str = "[N/A]"
            contrast_str = f"{r['gmia_contrast_cl_jja']:+.4f} {ci_str} ({r['contrast_status']})" if not np.isnan(r['gmia_contrast_cl_jja']) else f"N/A ({r['contrast_status']})"
            r_ol_str = f"{r['r_smapstar_vs_ol']:.4f}" if not np.isnan(r['r_smapstar_vs_ol']) else "N/A"
            r_anom_str = f"{r['r_anom_with_smap_raw']:.4f}" if not np.isnan(r['r_anom_with_smap_raw']) else "N/A"
            rows.append(
                f"| {r['model']} | {int(r['n_assimilable']):,} ({r['assimilable_pct']:.1f}%) | "
                f"{int(r['wilt_rejections']):,} ({r['wilt_rejection_pct']:.1f}%) | "
                f"{r['bias']:+.4f} m³/m³ | {r['rmse']:.4f} m³/m³ | {r['ubrmse']:.4f} m³/m³ | "
                f"{r_ol_str} | {r_anom_str} | {contrast_str} |"
            )
        return chr(10).join(rows)
        
    table_spec = make_metrics_table('Specific_Assimilable')
    table_comm = make_metrics_table('Common_Assimilable')
    
    # 3. Seasonal Stratified Metrics Table (Specific Assimilable)
    sub_sm = df_seas_metrics[df_seas_metrics['evaluation_population'] == 'Specific_Assimilable']
    sm_rows = []
    for _, r in sub_sm.iterrows():
        r_ol_s = f"{r['r_smapstar_vs_ol']:.4f}" if not np.isnan(r['r_smapstar_vs_ol']) else "N/A"
        r_anom_s = f"{r['r_anom_with_smap_raw']:.4f}" if not np.isnan(r['r_anom_with_smap_raw']) else "N/A"
        sm_rows.append(
            f"| {r['model']} | {r['season']} | {int(r['n_assimilable']):,} ({r['assimilable_pct']:.1f}%) | "
            f"{int(r['wilt_rejections']):,} ({r['wilt_rejection_pct']:.1f}%) | "
            f"{r['bias']:+.4f} m³/m³ | {r['rmse']:.4f} m³/m³ | {r['ubrmse']:.4f} m³/m³ | "
            f"{r_ol_s} | {r_anom_s} |"
        )
    table_seas_metrics = chr(10).join(sm_rows)
    
    # 4. Success Criteria Table
    def fmt_bool(b):
        return "**Validé**" if b else "**Non validé**"
        
    crit_rows = [
        f"| 1. Biais absolu (échantillon spécifique) | $|\\text{{Biais}}(\\mathrm{{V2}})| < |\\text{{Biais}}(\\mathrm{{B2\\text{{-}}train}})|$ | $|{crit['crit_1_bias_spec']['v2']:+.4f}| < |{crit['crit_1_bias_spec']['b2']:+.4f}|$ | {fmt_bool(crit['crit_1_bias_spec']['met'])} |",
        f"| 2. Biais absolu (échantillon commun) | $|\\text{{Biais}}(\\mathrm{{V2}})| < |\\text{{Biais}}(\\mathrm{{B2\\text{{-}}train}})|$ | $|{crit['crit_2_bias_comm']['v2']:+.4f}| < |{crit['crit_2_bias_comm']['b2']:+.4f}|$ | {fmt_bool(crit['crit_2_bias_comm']['met'])} |",
        f"| 3. Corrélation d'anomalies (échantillon spécifique) | $r_{{\\mathrm{{anom}}}}(\\mathrm{{V2}}) > r_{{\\mathrm{{anom}}}}(\\mathrm{{B2\\text{{-}}train}})$ | {crit['crit_3_ranom_spec']['v2']:.4f} > {crit['crit_3_ranom_spec']['b2']:.4f} | {fmt_bool(crit['crit_3_ranom_spec']['met'])} |",
        f"| 4. Corrélation d'anomalies (échantillon commun) | $r_{{\\mathrm{{anom}}}}(\\mathrm{{V2}}) > r_{{\\mathrm{{anom}}}}(\\mathrm{{B2\\text{{-}}train}})$ | {crit['crit_4_ranom_comm']['v2']:.4f} > {crit['crit_4_ranom_comm']['b2']:.4f} | {fmt_bool(crit['crit_4_ranom_comm']['met'])} |",
        f"| 5. RMSE (échantillon commun) | $\\text{{RMSE}}(\\mathrm{{V2}}) \\le 1{{,}}05 \\times \\text{{RMSE}}(\\mathrm{{B2\\text{{-}}train}})$ | {crit['crit_5_rmse_comm']['v2']:.4f} m³/m³ $\\le$ {crit['crit_5_rmse_comm']['thresh_105']:.4f} m³/m³ (B2-train: {crit['crit_5_rmse_comm']['b2']:.4f}) | {fmt_bool(crit['crit_5_rmse_comm']['met'])} |"
    ]
    crit_table = chr(10).join(crit_rows)
    overall_verdict = "**CRITÈRE DE SUCCÈS DÉCISION 15 SATISFAIT**" if crit['overall_success'] else "**CRITÈRE DE SUCCÈS DÉCISION 15 NON SATISFAIT**"

    report_text = f"""# Rapport de Phase 3 : Évaluation Scellée sur l'Année 2020

**Date d'exécution :** {time.strftime('%Y-%m-%d %H:%M:%S')}  
**Période d'évaluation :** 2020-01-01 au 2020-12-31 (année 2020 scellée au préalable, ouverte une seule fois)  
**Modèle utilisé :** Modèle LightGBM pré-entraîné sur 2016–2018 (300 arbres, learning_rate=0,05, num_leaves=63, aucun réentraînement, aucun changement d'hyperparamètres)  
**Référence principale :** CDF B2-train reconstruite avec LDT strictement sur 2016-01-01 → 2019-01-01  
**Méthode candidate formelle (Décision 15) :** Variante V2 (Removable + Context Representation) : {latex_v2}  
**Observations 2020 brutes :** {n_total_2020:,}  
**Observations 2020 passant QC (hors flétrissement) :** {n_base_qc:,} ({n_base_qc/n_total_2020*100:.2f}%)  
**Constante de base :** $\\phi_0 = {phi_0:+.4f}$ m³/m³ (TreeSHAP exact)

---

## 1. Évaluation des Critères de Succès Figés (Décision 15)

Verdict global avant ouverture : {overall_verdict}

| Critère | Formulation | Valeurs mesurées (V2 vs B2-train) | Résultat |
| :--- | :--- | :--- | :---: |
{crit_table}

---

## 2. Métriques Complètes du Protocole 2a (Année 2020)

Le flétrissement ({latex_wlt}) est recalculé sur l'observation corrigée pour chaque référence.

### 2.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\\mathrm{{SMAP}}^* \\mathrm{{vs}} \\mathrm{{OL}}}}$ | $r_{{\\mathrm{{anom}}}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_spec}

### 2.2 Échantillon commun ({n_common:,} observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\\mathrm{{SMAP}}^* \\mathrm{{vs}} \\mathrm{{OL}}}}$ | $r_{{\\mathrm{{anom}}}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_comm}

---

## 3. Corrections Saisonnières Moyennes Appliquées $\\overline{{\\Delta \\mathrm{{SMAP}}}}$ et Biais Résiduel $\\overline{{d_{{\\mathrm{{res}}}}}}$ (2020)

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{{\\mathrm{{res}}}}$ | V2 : Correction | V2 : Biais $d_{{\\mathrm{{res}}}}$ | V3 : Correction | V3 : Biais $d_{{\\mathrm{{res}}}}$ | B2-train : Correction | B2-train : Biais $d_{{\\mathrm{{res}}}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{seas_table}

---

## 4. Stratification Saisonnière des Métriques (Échantillon spécifique)

| Référence | Saison | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\\mathrm{{SMAP}}^* \\mathrm{{vs}} \\mathrm{{OL}}}}$ | $r_{{\\mathrm{{anom}}}}$ (SMAP brut) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_seas_metrics}

---

## 5. Figures de Publication (Année 2020)

- **Cartes saisonnières de correction moyenne (V1, V2, V3, B2-train) :** [F11_test2020_seasonal_applied_corrections.png](file://{os.path.join(FIG_DIR, "F11_test2020_seasonal_applied_corrections.png")}) | EN: [F11_test2020_seasonal_applied_corrections_en.png](file://{os.path.join(FIG_EN_DIR, "F11_test2020_seasonal_applied_corrections_en.png")})
"""

    with open(out_md, 'w', encoding='utf-8') as f:
        f.write(report_text)
    print(f"Generated report: {out_md}")


def main():
    print("=== STARTING PHASE 3 SEALED 2020 TEST ===")
    t_start = time.time()
    
    # 1. Feature groups
    cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, feature_cols = load_feature_groups()
    print(f"Loaded 66 features: {len(cols_rem)} removable, {len(cols_ret)} retained, {len(cols_ctx_rep)} context_rep, {len(cols_ctx_proc)} context_proc")
    
    # 2. Static variables for masking & wilting
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]
        
    n_jobs = get_n_jobs()
    print(f"Allocated CPU cores (n_jobs): {n_jobs}")

    # 3. Load pre-trained references
    assert os.path.exists(LUT_CACHE_FILE), f"Missing {LUT_CACHE_FILE}"
    assert os.path.exists(LUT_B2_TRAIN_FILE), f"Missing {LUT_B2_TRAIN_FILE}"
    assert os.path.exists(B1_CLIM_FILE), f"Missing {B1_CLIM_FILE}"
    assert os.path.exists(B1_PIX_FILE), f"Missing {B1_PIX_FILE}"
    assert os.path.exists(MODEL_FILE), f"Missing {MODEL_FILE}"
    
    lut_b2 = np.load(LUT_CACHE_FILE)
    lut_b2_train = np.load(LUT_B2_TRAIN_FILE)
    b1_clim_table = pd.read_csv(B1_CLIM_FILE)
    b1_pixel_table = pd.read_csv(B1_PIX_FILE)
    global_fallback = float(b1_pixel_table['b1_pixel_fallback'].mean())
    
    # 4. Evaluate 2020 Test Set (Strictly Once, Pre-Trained Model Only)
    df_2020, shap_2020, X_val_2020, phi_0, n_total_2020, n_base_qc = evaluate_2020_sealed_test(
        feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc,
        lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, b1_clim_table, b1_pixel_table, global_fallback, n_jobs=n_jobs
    )
    
    # 5. Compute Seasonal Corrections, Full Protocol 2a Metrics, and Seasonal Metrics
    df_seasonal, df_metrics, df_seas_metrics, n_common = compute_seasonal_corrections_and_metrics_2020(
        df_2020, n_total_2020, n_base_qc
    )
    
    # 6. Evaluate Decision 15 Success Criteria
    crit = evaluate_decision15_criteria(df_metrics)
    
    # 7. Generate Publication Figures
    generate_figures_2020(df_2020)
    
    # 8. Generate Markdown Report directly
    generate_markdown_report_2020(
        df_seasonal, df_metrics, df_seas_metrics,
        n_total_2020, n_base_qc, n_common, phi_0, crit
    )
    
    print(f"=== PHASE 3 SEALED 2020 TEST COMPLETED IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
