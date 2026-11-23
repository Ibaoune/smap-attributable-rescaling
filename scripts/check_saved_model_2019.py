#!/usr/bin/env python3
"""
Contrôle du modèle sauvegardé data/interim/lgb_model_2016_2018.txt sur l'année 2019.

Ce script :
1. Charge le modèle sauvegardé data/interim/lgb_model_2016_2018.txt sans aucun réentraînement.
2. Évalue le modèle sur l'année de validation 2019 avec l'évaluation exacte du script 10
   (Lundberg TreeSHAP C++ API pred_contrib=True, regroupements V1/V2/V3, baselines B0, B1, B2, B2-train).
3. Compare chaque métrique aux valeurs publiées dans reports/phase2c_lgb_metrics.csv.
4. Écrit l'écart par métrique dans reports/check_saved_model_2019.csv.
"""

import os
import sys
import time
import yaml
import netCDF4 as nc
import numpy as np
import pandas as pd
import lightgbm as lgb

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
INTERIM_DIR = os.path.join(BASE_DIR, "data", "interim")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
CONFIG_PREDICTORS = os.path.join(BASE_DIR, "configs", "predictors.yaml")

STATIC_FILE = os.path.join(INTERIM_DIR, "static", "static_predictors.nc")
CDF_OBS_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "smap_obs_cdf_NorthMor.nc")
CDF_B2_TRAIN_OBS = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_obs", "smap_obs_cdf_NorthMor.nc")
LUT_CACHE_FILE = os.path.join(INTERIM_DIR, "lut_cdf_rescaling.npy")
LUT_B2_TRAIN_FILE = os.path.join(INTERIM_DIR, "lut_cdf_b2_train.npy")
B1_CLIM_FILE = os.path.join(INTERIM_DIR, "b1_climatology_2016_2018.csv")
B1_PIX_FILE = os.path.join(INTERIM_DIR, "b1_pixel_fallback_2016_2018.csv")
MODEL_FILE = os.path.join(INTERIM_DIR, "lgb_model_2016_2018.txt")
PHASE2C_METRICS_FILE = os.path.join(REPORTS_DIR, "phase2c_lgb_metrics.csv")
CHECK_OUTPUT_FILE = os.path.join(REPORTS_DIR, "check_saved_model_2019.csv")

TEXTURE_NAMES = {
    1: 'Sand', 2: 'Loamy Sand', 3: 'Sandy Loam', 4: 'Silt Loam', 5: 'Silt',
    6: 'Loam', 7: 'Sandy Clay Loam', 8: 'Silty Clay Loam', 9: 'Clay Loam',
    10: 'Sandy Clay', 11: 'Silty Clay', 12: 'Clay', 13: 'Organic Material',
    14: 'Water', 15: 'Bedrock', 16: 'Other'
}


def get_n_jobs():
    if "SLURM_CPUS_PER_TASK" in os.environ:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    return os.cpu_count() or 8


def load_feature_groups():
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


def main():
    print("=== STARTING CHECK OF SAVED LIGHTGBM MODEL ON 2019 VALIDATION SET ===")
    t0 = time.time()
    n_jobs = get_n_jobs()
    print(f"Allocated CPUs: {n_jobs}")

    assert os.path.exists(MODEL_FILE), f"Missing model file: {MODEL_FILE}"
    assert os.path.exists(PHASE2C_METRICS_FILE), f"Missing reference metrics file: {PHASE2C_METRICS_FILE}"

    # 1. Feature groups
    cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, feature_cols = load_feature_groups()
    idx_rem = [feature_cols.index(c) for c in cols_rem]
    idx_ret = [feature_cols.index(c) for c in cols_ret]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_ctx_rep]
    idx_ctx_proc = [feature_cols.index(c) for c in cols_ctx_proc]

    # 2. Static masks
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]

    # 3. Reference LUTs & B1 Climatology
    lut_b2 = np.load(LUT_CACHE_FILE)
    lut_b2_train = np.load(LUT_B2_TRAIN_FILE)
    b1_clim_table = pd.read_csv(B1_CLIM_FILE)
    b1_pixel_table = pd.read_csv(B1_PIX_FILE)

    # Script 10 global fallback definition (mean of target_d across 2016-2018 non-wilt QC)
    qc_cols = [
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    dfs = []
    for y in [2016, 2017, 2018]:
        pq = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(pq, columns=['target_d'] + qc_cols)
        qc_nw = (
            (df_y['qc_lis_rain'] == 0) &
            (df_y['qc_lis_frozen'] == 0) &
            (df_y['qc_lis_forest'] == 0) &
            (df_y['qc_lis_snow'] == 0) &
            (df_y['qc_lis_high_veg'] == 0) &
            (df_y['qc_lis_sat_limit'] == 0) &
            (df_y['qc_lis_vwc_gt_5'] == 0) &
            (~df_y['target_d'].isna())
        )
        dfs.append(df_y.loc[qc_nw, ['target_d']])
    full_train = pd.concat(dfs, ignore_index=True)
    global_fallback = float(full_train['target_d'].mean())
    print(f"Script 10 global_fallback calculated: {global_fallback:.6f}")

    # 4. Load 2019 validation set
    pq_2019 = os.path.join(TABLE_DIR, "sar_modeling_table_2019.parquet")
    meta_cols = [
        'pixel_id', 'month', 'season', 'gmia_class', 'spatial_block_id',
        'smap_sm', 'ssm_ol_daily', 'target_d', 'split',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    load_cols = list(set(feature_cols + meta_cols))
    df_raw_2019 = pd.read_parquet(pq_2019, columns=load_cols)
    n_total_2019 = len(df_raw_2019)

    qc_nw = (
        (df_raw_2019['qc_lis_rain'] == 0) &
        (df_raw_2019['qc_lis_frozen'] == 0) &
        (df_raw_2019['qc_lis_forest'] == 0) &
        (df_raw_2019['qc_lis_snow'] == 0) &
        (df_raw_2019['qc_lis_high_veg'] == 0) &
        (df_raw_2019['qc_lis_sat_limit'] == 0) &
        (df_raw_2019['qc_lis_vwc_gt_5'] == 0) &
        (~df_raw_2019['smap_sm'].isna())
    )
    df_2019 = df_raw_2019[qc_nw].copy().reset_index(drop=True)
    n_base_qc = len(df_2019)
    print(f"Loaded 2019 non-wilt QC: {n_base_qc:,} rows")

    r_idx = df_2019['pixel_id'].values // 181
    c_idx = df_2019['pixel_id'].values % 181
    wlt_thresh = soil_wlt[r_idx, c_idx] + 0.02
    soil_ids = dom_soil[r_idx, c_idx]
    df_2019['texture_name'] = [TEXTURE_NAMES.get(int(s), "Other") for s in soil_ids]
    df_2019['wlt_thresh'] = wlt_thresh

    tile_map = np.full((121, 181), -1, dtype=int)
    t = 0
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                tile_map[r, c] = t
                t += 1
    tiles = tile_map[r_idx, c_idx]

    # 5. Baseline B0
    df_2019['smap_star_b0'] = df_2019['smap_sm']
    df_2019['d_res_b0'] = df_2019['target_d']
    df_2019['under_wlt_b0'] = (df_2019['smap_star_b0'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b0'] = (
        (df_2019['under_wlt_b0'] == 0) &
        (df_2019['smap_star_b0'] >= 0.01) &
        (df_2019['smap_star_b0'] <= 1.0) &
        (~df_2019['smap_star_b0'].isna())
    )

    # 6. Baseline B1
    df_2019 = df_2019.merge(b1_clim_table, on=['pixel_id', 'month'], how='left')
    df_2019 = df_2019.merge(b1_pixel_table, on='pixel_id', how='left')
    c_hat = df_2019['b1_bias'].fillna(df_2019['b1_pixel_fallback']).fillna(global_fallback).values.astype(np.float32)
    df_2019.drop(columns=['b1_bias', 'b1_pixel_fallback'], inplace=True)

    df_2019['smap_star_b1'] = df_2019['smap_sm'] - c_hat
    df_2019['d_res_b1'] = df_2019['target_d'] - c_hat
    df_2019['under_wlt_b1'] = (df_2019['smap_star_b1'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b1'] = (
        (df_2019['under_wlt_b1'] == 0) &
        (df_2019['smap_star_b1'] >= 0.01) &
        (df_2019['smap_star_b1'] <= 1.0) &
        (~df_2019['smap_star_b1'].isna())
    )

    # 7. Baseline B2 (Legacy) & B2-train
    months = df_2019['month'].values
    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta = obs_x[1, :, :] - obs_x[0, :, :]
    obs_x0 = obs_x[0, :, :]
    m_indices = months - 1
    smap_vals = df_2019['smap_sm'].values

    def apply_lut(lut_arr, obs_x0_arr, obs_delta_arr):
        res = np.full(len(df_2019), np.nan, dtype=np.float32)
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
    df_2019['smap_star_b2'] = smap_star_b2
    df_2019['b2_has_value'] = (~np.isnan(smap_star_b2)).astype(int)
    df_2019['d_res_b2'] = smap_star_b2 - df_2019['ssm_ol_daily']
    df_2019['under_wlt_b2'] = (df_2019['smap_star_b2'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b2'] = (
        (df_2019['b2_has_value'] == 1) &
        (df_2019['under_wlt_b2'] == 0) &
        (df_2019['smap_star_b2'] >= 0.01) &
        (df_2019['smap_star_b2'] <= 1.0)
    )

    with nc.Dataset(CDF_B2_TRAIN_OBS) as ds_o_tr:
        obs_x_tr = ds_o_tr.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta_tr = obs_x_tr[1, :, :] - obs_x_tr[0, :, :]
    obs_x0_tr = obs_x_tr[0, :, :]

    smap_star_b2_tr = apply_lut(lut_b2_train, obs_x0_tr, obs_delta_tr)
    df_2019['smap_star_b2_train'] = smap_star_b2_tr
    df_2019['b2_train_has_value'] = (~np.isnan(smap_star_b2_tr)).astype(int)
    df_2019['d_res_b2_train'] = smap_star_b2_tr - df_2019['ssm_ol_daily']
    df_2019['under_wlt_b2_train'] = (df_2019['smap_star_b2_train'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b2_train'] = (
        (df_2019['b2_train_has_value'] == 1) &
        (df_2019['under_wlt_b2_train'] == 0) &
        (df_2019['smap_star_b2_train'] >= 0.01) &
        (df_2019['smap_star_b2_train'] <= 1.0)
    )

    # 8. Apply saved LightGBM model WITHOUT retraining
    print(f"Loading pre-trained LightGBM model from {MODEL_FILE}...")
    booster = lgb.Booster(model_file=MODEL_FILE)
    print(f"Booster loaded successfully: {booster.num_trees()} trees, {booster.num_feature()} features")

    X_val_2019 = df_2019[feature_cols].values.astype(np.float32)
    print("Computing exact Lundberg TreeSHAP via Booster.predict(..., pred_contrib=True)...")
    t_shap = time.time()
    contrib_2019 = booster.predict(X_val_2019, pred_contrib=True, num_threads=n_jobs)
    print(f"Exact TreeSHAP computed in {time.time()-t_shap:.1f}s")

    phi_0 = float(contrib_2019[0, -1])
    shap_2019 = contrib_2019[:, :-1]

    # Group sums
    phi_rem = shap_2019[:, idx_rem].sum(axis=1)
    phi_ret = shap_2019[:, idx_ret].sum(axis=1)
    phi_ctx_rep = shap_2019[:, idx_ctx_rep].sum(axis=1)
    phi_ctx_proc = shap_2019[:, idx_ctx_proc].sum(axis=1)

    # V1 = removable
    df_2019['smap_star_v1'] = df_2019['smap_sm'] - (phi_rem + phi_0)
    df_2019['d_res_v1'] = df_2019['smap_star_v1'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v1'] = (df_2019['smap_star_v1'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v1'] = (
        (df_2019['under_wlt_v1'] == 0) &
        (df_2019['smap_star_v1'] >= 0.01) &
        (df_2019['smap_star_v1'] <= 1.0) &
        (~df_2019['smap_star_v1'].isna())
    )

    # V2 = removable + context_rep
    df_2019['smap_star_v2'] = df_2019['smap_sm'] - (phi_rem + phi_ctx_rep + phi_0)
    df_2019['d_res_v2'] = df_2019['smap_star_v2'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v2'] = (df_2019['smap_star_v2'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v2'] = (
        (df_2019['under_wlt_v2'] == 0) &
        (df_2019['smap_star_v2'] >= 0.01) &
        (df_2019['smap_star_v2'] <= 1.0) &
        (~df_2019['smap_star_v2'].isna())
    )

    # V3 = removable + all context
    df_2019['smap_star_v3'] = df_2019['smap_sm'] - (phi_rem + phi_ctx_rep + phi_ctx_proc + phi_0)
    df_2019['d_res_v3'] = df_2019['smap_star_v3'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v3'] = (df_2019['smap_star_v3'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v3'] = (
        (df_2019['under_wlt_v3'] == 0) &
        (df_2019['smap_star_v3'] >= 0.01) &
        (df_2019['smap_star_v3'] <= 1.0) &
        (~df_2019['smap_star_v3'].isna())
    )

    # Common sample
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b0'] &
        df_2019['assimilable_b1'] &
        df_2019['assimilable_b2'] &
        df_2019['assimilable_b2_train'] &
        df_2019['assimilable_v1'] &
        df_2019['assimilable_v2'] &
        df_2019['assimilable_v3']
    )

    # 9. Compute metrics
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
                sub = df_2019[df_2019[assim_col]].copy()
                n_wlt = int(df_2019[wlt_col].sum())
            else:
                sub = df_2019[df_2019[pop_col]].copy()
                n_wlt = 0
                
            n_assim = len(sub)
            d_vals = sub[d_col].values
            bias = float(np.mean(d_vals))
            rmse = float(np.sqrt(np.mean(d_vals ** 2)))
            ubrmse = float(np.sqrt(np.var(d_vals)))
            
            ol_vals = sub['ssm_ol_daily'].values
            star_vals = sub[star_col].values
            r_ol = float(np.corrcoef(star_vals, ol_vals)[0, 1]) if len(sub) > 1 else np.nan
            
            anom_raw = (sub['smap_sm'] - sub.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')).values
            anom_star = (sub[star_col] - sub.groupby(['pixel_id', 'month'])[star_col].transform('mean')).values
            valid_a = ~(np.isnan(anom_raw) | np.isnan(anom_star))
            r_anom = float(np.corrcoef(anom_raw[valid_a], anom_star[valid_a])[0, 1]) if valid_a.sum() > 10 else np.nan
            
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

    df_check = pd.DataFrame(metrics_rows)

    # 10. Compare with phase2c_lgb_metrics.csv
    print(f"Loading reference metrics from {PHASE2C_METRICS_FILE}...")
    df_ref = pd.read_csv(PHASE2C_METRICS_FILE)

    numeric_cols = [
        'n_assimilable', 'assimilable_pct', 'wilt_rejections', 'wilt_rejection_pct',
        'bias', 'rmse', 'ubrmse', 'r_smapstar_vs_ol', 'r_anom_with_smap_raw',
        'gmia_contrast_cl_jja', 'gmia_contrast_ci_low', 'gmia_contrast_ci_high',
        'n_px_gt30', 'n_px_0', 'n_blocks_gt30', 'n_blocks_0'
    ]

    diff_rows = []
    max_abs_diff_by_model = {}

    for pop in ['Specific_Assimilable', 'Common_Assimilable']:
        for name, _, _, _, _ in models:
            sub_ref = df_ref[(df_ref['evaluation_population'] == pop) & (df_ref['model'] == name)].iloc[0]
            sub_chk = df_check[(df_check['evaluation_population'] == pop) & (df_check['model'] == name)].iloc[0]
            
            for m in numeric_cols:
                v_ref = float(sub_ref[m]) if pd.notnull(sub_ref[m]) else np.nan
                v_chk = float(sub_chk[m]) if pd.notnull(sub_chk[m]) else np.nan
                diff = v_chk - v_ref if (pd.notnull(v_ref) and pd.notnull(v_chk)) else (0.0 if pd.isnull(v_ref) and pd.isnull(v_chk) else np.nan)
                
                diff_rows.append({
                    'evaluation_population': pop,
                    'model': name,
                    'metric': m,
                    'phase2c_val': v_ref,
                    'saved_model_val': v_chk,
                    'diff': diff
                })
                
                if pd.notnull(diff):
                    k = f"{pop}_{name}"
                    max_abs_diff_by_model[k] = max(max_abs_diff_by_model.get(k, 0.0), abs(diff))

    df_diffs = pd.DataFrame(diff_rows)
    df_diffs.to_csv(CHECK_OUTPUT_FILE, index=False)
    print(f"Saved check discrepancies to: {CHECK_OUTPUT_FILE}")

    print("\n--- Summary of Max Absolute Discrepancies by Population & Model ---")
    for k, v in max_abs_diff_by_model.items():
        print(f"  {k:45s}: max |diff| = {v:.2e}")

    # Specific check on candidate model LGB_V2_Rem_ContextRep
    v2_diffs = df_diffs[df_diffs['model'] == 'LGB_V2_Rem_ContextRep']
    max_v2_diff = float(v2_diffs['diff'].abs().max())
    print(f"\nCandidate Model LGB_V2_Rem_ContextRep max absolute discrepancy: {max_v2_diff:.2e}")
    if max_v2_diff < 1e-12:
        print("PERFECT MATCH: The saved LightGBM model identically reproduces phase2c metrics to machine precision!")
    else:
        print("WARNING: Discrepancies observed!")

    print(f"\nTotal check execution time: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
