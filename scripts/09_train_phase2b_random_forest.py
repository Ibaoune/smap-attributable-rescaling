#!/usr/bin/env python3
"""
scripts/09_train_phase2b_random_forest.py

PHASE 2b : First Machine Learning Model - Random Forest (RF).
Strict protocol:
- Training: 2016-2018 (split == 'train'), non-wilt QC valid observations.
- Stratified sampling: 20 observations per pixel across 2016-2018 (uniform spatial representation).
- Evaluation: 2019 (split == 'val').
- 2020 remains strictly sealed!
- Predictor groups from configs/predictors.yaml (removable, retained, context).
- 5-fold balanced spatial cross-validation + TreeSHAP stability (Jalon H4).
- SAR Attribution correction: SMAP* = SMAP - (Phi_removable + phi_0).
  Sensitivity variant: SMAP* = SMAP - (Phi_removable + Phi_context + phi_0).
- Recalculated wilt rule (SMAP* < theta_wlt + 0.02) and identical protocol 2a evaluation.
- Zero interpretive phrasing (only computed facts).
"""

import os
import time
import yaml
import netCDF4 as nc
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import shap

# Paths (Relative to repository root)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
CDF_MODEL_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "noahmp_model_cdf_NorthMor.nc")
CDF_OBS_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "smap_obs_cdf_NorthMor.nc")
LUT_CACHE_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_rescaling.npy")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def get_feature_groups():
    """Define feature column lists for the three groups."""
    cols_removable = [
        'elevation', 'slope_deg', 'slopetype',
        'vegetation_opacity', 'vegetation_water_content',
        'surface_temperature', 'roughness_coefficient', 'roughness_is_nan',
        'static_water_body_fraction', 'surface_flag', 'pass_is_pm'
    ] + [f'texture_frac_{i}' for i in range(1, 17)] + [f'landcover_frac_{i}' for i in range(1, 21)]

    cols_retained = [
        'irrigfrac_gmia', 'delta_precip_imerg_mswep', 'precip_mswep',
        'precip_1d_sum', 'precip_3d_sum', 'precip_7d_sum', 'days_since_rain_gt_1mm'
    ]

    cols_context = [
        'doy_sin', 'doy_cos', 'ssm_lag1', 'sm_layer2_lag1', 'sm_layer3_lag1', 'sm_layer4_lag1',
        'rzsm_lag1', 'et_lag1', 'surface_runoff_lag1', 'baseflow_lag1', 'lai_lag1', 'ssm_rel_wilt_ref_lag1'
    ]
    return cols_removable, cols_retained, cols_context


def load_training_data(feature_cols):
    """Load 2016-2018 non-wilt QC valid data and perform stratified sampling per pixel."""
    print("--- Loading 2016-2018 Training Data (Stratified Sampling per Pixel) ---")
    t0 = time.time()
    
    meta_cols = ['pixel_id', 'date', 'month', 'season', 'pass', 'spatial_block_id', 'spatial_fold', 'split', 'target_d', 'smap_sm', 'ssm_ol_daily']
    qc_cols = ['qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow', 'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5']
    load_cols = list(set(feature_cols + meta_cols + qc_cols))
    
    dfs = []
    for y in [2016, 2017, 2018]:
        p = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df = pd.read_parquet(p, columns=load_cols)
        qc_nw = (
            (df['split'] == 'train') &
            (df['qc_lis_rain'] == 0) &
            (df['qc_lis_frozen'] == 0) &
            (df['qc_lis_forest'] == 0) &
            (df['qc_lis_snow'] == 0) &
            (df['qc_lis_high_veg'] == 0) &
            (df['qc_lis_sat_limit'] == 0) &
            (df['qc_lis_vwc_gt_5'] == 0) &
            (~df['smap_sm'].isna())
        )
        df_valid = df[qc_nw].copy()
        df_valid['rand'] = np.random.RandomState(42 + y).rand(len(df_valid))
        sampled = df_valid.sort_values(['pixel_id', 'rand']).groupby('pixel_id').head(10)
        dfs.append(sampled)
        print(f"  Year {y}: loaded & sampled {len(sampled):,} rows")
        
    train_df = pd.concat(dfs, ignore_index=True)
    train_df['rand'] = np.random.RandomState(42).rand(len(train_df))
    train_sampled = train_df.sort_values(['pixel_id', 'rand']).groupby('pixel_id').head(20).reset_index(drop=True)
    
    # Impute roughness_coefficient NaN with median/0.0
    train_sampled['roughness_coefficient'] = train_sampled['roughness_coefficient'].fillna(0.0)
    
    n_train = len(train_sampled)
    n_pixels = int(train_sampled['pixel_id'].nunique())
    print(f"Training set ready: {n_train:,} rows across {n_pixels:,} pixels in {time.time()-t0:.2f}s")
    return train_sampled


def run_spatial_cross_validation(train_df, feature_cols, cols_removable, cols_retained, cols_context):
    """Perform 5-fold balanced spatial cross-validation and evaluate SHAP stability (Jalon H4)."""
    print("--- Running 5-Fold Spatial Cross-Validation and SHAP Stability (H4) ---")
    t0 = time.time()
    
    idx_rem = [feature_cols.index(c) for c in cols_removable]
    idx_ret = [feature_cols.index(c) for c in cols_retained]
    idx_ctx = [feature_cols.index(c) for c in cols_context]
    
    X_train_full = train_df[feature_cols].values.astype(np.float32)
    y_train_full = train_df['target_d'].values.astype(np.float32)
    folds = train_df['spatial_fold'].values
    
    cv_records = []
    shap_stability_records = []
    
    for k in range(5):
        t_fold = time.time()
        train_mask = (folds != k)
        val_mask = (folds == k)
        
        X_tr, y_tr = X_train_full[train_mask], y_train_full[train_mask]
        X_va, y_va = X_train_full[val_mask], y_train_full[val_mask]
        
        rf = RandomForestRegressor(
            n_estimators=50, max_depth=14, min_samples_leaf=10, max_features='sqrt',
            n_jobs=8, random_state=42 + k
        )
        rf.fit(X_tr, y_tr)
        
        # OOF Predictions
        y_pred = rf.predict(X_va)
        rmse_oof = float(np.sqrt(np.mean((y_va - y_pred) ** 2)))
        mae_oof = float(np.mean(np.abs(y_va - y_pred)))
        bias_oof = float(np.mean(y_pred - y_va))
        r2_oof = float(1.0 - np.sum((y_va - y_pred) ** 2) / np.sum((y_va - np.mean(y_va)) ** 2))
        
        # SHAP on validation fold sample (exact TreeSHAP on 500 samples)
        np.random.seed(42 + k)
        sample_idx = np.random.choice(len(X_va), size=min(500, len(X_va)), replace=False)
        X_sample = X_va[sample_idx]
        
        explainer = shap.TreeExplainer(rf)
        sv = explainer.shap_values(X_sample, check_additivity=False)
        
        # Mean absolute attribution per group
        mean_abs_rem = float(np.mean(np.abs(sv[:, idx_rem].sum(axis=1))))
        mean_abs_ret = float(np.mean(np.abs(sv[:, idx_ret].sum(axis=1))))
        mean_abs_ctx = float(np.mean(np.abs(sv[:, idx_ctx].sum(axis=1))))
        
        tot_attr = mean_abs_rem + mean_abs_ret + mean_abs_ctx
        pct_rem = mean_abs_rem / tot_attr * 100.0 if tot_attr > 0 else 0.0
        pct_ret = mean_abs_ret / tot_attr * 100.0 if tot_attr > 0 else 0.0
        pct_ctx = mean_abs_ctx / tot_attr * 100.0 if tot_attr > 0 else 0.0
        
        cv_records.append({
            'spatial_fold': k,
            'n_train': int(train_mask.sum()),
            'n_val': int(val_mask.sum()),
            'oof_rmse': rmse_oof,
            'oof_mae': mae_oof,
            'oof_bias': bias_oof,
            'oof_r2': r2_oof,
            'fit_time_s': time.time() - t_fold
        })
        
        shap_stability_records.append({
            'spatial_fold': k,
            'mean_abs_removable': mean_abs_rem,
            'mean_abs_retained': mean_abs_ret,
            'mean_abs_context': mean_abs_ctx,
            'share_removable_pct': pct_rem,
            'share_retained_pct': pct_ret,
            'share_context_pct': pct_ctx
        })
        print(f"  Fold {k}: OOF RMSE = {rmse_oof:.4f}, R2 = {r2_oof:.4f} | Attributions: Rem={pct_rem:.1f}%, Ret={pct_ret:.1f}%, Ctx={pct_ctx:.1f}% ({time.time()-t_fold:.1f}s)")
        
    df_cv = pd.DataFrame(cv_records)
    df_shap_stab = pd.DataFrame(shap_stability_records)
    print(f"5-Fold CV completed in {time.time()-t0:.1f}s | Mean OOF RMSE: {df_cv['oof_rmse'].mean():.4f}, Mean R2: {df_cv['oof_r2'].mean():.4f}")
    return df_cv, df_shap_stab


def train_full_model(train_df, feature_cols):
    """Train final Random Forest on full 2016-2018 training sample."""
    print("--- Training Full Random Forest Model on 2016-2018 ---")
    t0 = time.time()
    X_train = train_df[feature_cols].values.astype(np.float32)
    y_train = train_df['target_d'].values.astype(np.float32)
    
    rf_full = RandomForestRegressor(
        n_estimators=50, max_depth=14, min_samples_leaf=10, max_features='sqrt',
        n_jobs=8, random_state=42
    )
    rf_full.fit(X_train, y_train)
    
    explainer = shap.TreeExplainer(rf_full)
    base_val = float(np.asarray(explainer.expected_value).flatten()[0])
    print(f"Full model trained in {time.time()-t0:.2f}s | Base value (phi_0): {base_val:+.6f}")
    return rf_full, explainer, base_val


def load_lut_and_climatology(landmask):
    """Load cached LUT and learn B1 climatology on 2016-2018 non-wilt QC."""
    lut = np.load(LUT_CACHE_FILE)
    
    dfs = []
    qc_cols = [
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    for y in [2016, 2017, 2018]:
        pq = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(pq, columns=['split', 'pixel_id', 'month', 'target_d'] + qc_cols)
        qc_nw = (
            (df_y['split'] == 'train') &
            (df_y['qc_lis_rain'] == 0) &
            (df_y['qc_lis_frozen'] == 0) &
            (df_y['qc_lis_forest'] == 0) &
            (df_y['qc_lis_snow'] == 0) &
            (df_y['qc_lis_high_veg'] == 0) &
            (df_y['qc_lis_sat_limit'] == 0) &
            (df_y['qc_lis_vwc_gt_5'] == 0)
        )
        dfs.append(df_y.loc[qc_nw, ['pixel_id', 'month', 'target_d']])
    train_df = pd.concat(dfs, ignore_index=True)
    climatology = train_df.groupby(['pixel_id', 'month'])['target_d'].mean().to_dict()
    pixel_fallback = train_df.groupby('pixel_id')['target_d'].mean().to_dict()
    global_fallback = float(train_df['target_d'].mean())
    return lut, climatology, pixel_fallback, global_fallback


def evaluate_phase2b_on_2019(rf_full, explainer, base_val, feature_cols, cols_removable, cols_retained, cols_context,
                             lut, climatology, pixel_fallback, global_fallback, landmask, soil_wlt, dom_soil):
    """Evaluate RF (SAR main and context sensitivity) alongside B0, B1, B2 on 2019 with EXACT protocol 2a."""
    print("--- Evaluating Phase 2b Models & Baselines on 2019 (Val Split) ---")
    t0 = time.time()
    
    idx_rem = [feature_cols.index(c) for c in cols_removable]
    idx_ret = [feature_cols.index(c) for c in cols_retained]
    idx_ctx = [feature_cols.index(c) for c in cols_context]
    
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
    
    # Non-wilt QC population
    qc_nw = (
        (df_raw_2019['qc_lis_rain'] == 0) &
        (df_raw_2019['qc_lis_frozen'] == 0) &
        (df_raw_2019['qc_lis_forest'] == 0) &
        (df_raw_2019['qc_lis_snow'] == 0) &
        (df_raw_2019['qc_lis_high_veg'] == 0) &
        (df_raw_2019['qc_lis_sat_limit'] == 0) &
        (df_raw_2019['qc_lis_vwc_gt_5'] == 0)
    )
    df_2019 = df_raw_2019[qc_nw].copy().reset_index(drop=True)
    n_base_qc = len(df_2019)
    print(f"Loaded 2019 non-wilt QC: {n_base_qc:,} rows ({n_base_qc/n_total_2019*100:.2f}%) in {time.time()-t0:.1f}s")
    
    # Impute roughness_coefficient
    df_2019['roughness_coefficient'] = df_2019['roughness_coefficient'].fillna(0.0)
    
    # Extract metadata arrays
    r_idx = df_2019['pixel_id'].values // 181
    c_idx = df_2019['pixel_id'].values % 181
    wlt_thresh = soil_wlt[r_idx, c_idx] + 0.02
    soil_ids = dom_soil[r_idx, c_idx]
    df_2019['texture_name'] = [TEXTURE_NAMES.get(s, "Other") for s in soil_ids]
    df_2019['wlt_thresh'] = wlt_thresh
    
    # Tile map for CDF
    tile_map = np.full((121, 181), -1, dtype=int)
    t = 0
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                tile_map[r, c] = t
                t += 1
    tiles = tile_map[r_idx, c_idx]
    
    # --- 1. Baseline B0: No correction ---
    df_2019['smap_star_b0'] = df_2019['smap_sm']
    df_2019['d_res_b0'] = df_2019['target_d']
    df_2019['under_wlt_b0'] = (df_2019['smap_star_b0'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b0'] = (
        (df_2019['under_wlt_b0'] == 0) &
        (df_2019['smap_star_b0'] >= 0.01) &
        (df_2019['smap_star_b0'] <= 1.0) &
        (~df_2019['smap_star_b0'].isna())
    )
    
    # --- 2. Baseline B1: Climatology ---
    p_ids = df_2019['pixel_id'].values
    months = df_2019['month'].values
    c_hat = np.array([climatology.get((p, m), pixel_fallback.get(p, global_fallback)) for p, m in zip(p_ids, months)], dtype=np.float32)
    df_2019['smap_star_b1'] = df_2019['smap_sm'] - c_hat
    df_2019['d_res_b1'] = df_2019['target_d'] - c_hat
    df_2019['under_wlt_b1'] = (df_2019['smap_star_b1'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b1'] = (
        (df_2019['under_wlt_b1'] == 0) &
        (df_2019['smap_star_b1'] >= 0.01) &
        (df_2019['smap_star_b1'] <= 1.0) &
        (~df_2019['smap_star_b1'].isna())
    )
    
    # --- 3. Baseline B2: CDF v2 offline ---
    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta = obs_x[1, :, :] - obs_x[0, :, :]
    obs_x0 = obs_x[0, :, :]
    m_indices = months - 1
    smap_vals = df_2019['smap_sm'].values
    smap_star_b2 = np.full(len(df_2019), np.nan, dtype=np.float32)
    
    for m in range(12):
        mask_m = (m_indices == m)
        if not np.any(mask_m):
            continue
        sub_tiles = tiles[mask_m]
        sub_smap = smap_vals[mask_m]
        sub_delta = obs_delta[m, sub_tiles]
        sub_x0 = obs_x0[m, sub_tiles]
        valid_sub = (sub_delta > 0) & (sub_tiles >= 0) & (sub_smap >= 0.01) & (sub_smap <= 1.0)
        binvals = np.zeros(len(sub_tiles), dtype=int)
        binvals[valid_sub] = np.floor((sub_smap[valid_sub] - sub_x0[valid_sub]) / sub_delta[valid_sub] + 0.5).astype(int) + 1
        binvals = np.clip(binvals, 1, 100) - 1
        rescaled_m = np.full(len(sub_tiles), np.nan, dtype=np.float32)
        rescaled_m[valid_sub] = lut[binvals[valid_sub], m, sub_tiles[valid_sub]]
        smap_star_b2[mask_m] = rescaled_m
        
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
    
    # --- 4. Model RF (SAR Main) and RF-ctx (Sensitivity Variant) ---
    print("Computing Random Forest TreeSHAP attributions on 2019 validation set...")
    t_shap0 = time.time()
    X_val_2019 = df_2019[feature_cols].values.astype(np.float32)
    
    # Compute additive tree path attributions
    sv_2019 = explainer.shap_values(X_val_2019, approximate=True, check_additivity=False)
    print(f"Computed TreeSHAP attributions on {len(df_2019):,} rows in {time.time()-t_shap0:.1f}s")
    
    phi_rem = sv_2019[:, idx_rem].sum(axis=1)
    phi_ret = sv_2019[:, idx_ret].sum(axis=1)
    phi_ctx = sv_2019[:, idx_ctx].sum(axis=1)
    
    df_2019['phi_removable'] = phi_rem
    df_2019['phi_retained'] = phi_ret
    df_2019['phi_context'] = phi_ctx
    
    # RF Main: correction = SMAP - (phi_rem + base_val)
    # The retained group (GMIA, precipitation discrepancies) is RETAINED in the observation!
    df_2019['smap_star_rf'] = df_2019['smap_sm'] - (phi_rem + base_val)
    df_2019['d_res_rf'] = df_2019['smap_star_rf'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_rf'] = (df_2019['smap_star_rf'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_rf'] = (
        (df_2019['under_wlt_rf'] == 0) &
        (df_2019['smap_star_rf'] >= 0.01) &
        (df_2019['smap_star_rf'] <= 1.0) &
        (~df_2019['smap_star_rf'].isna())
    )
    
    # RF-ctx (Sensitivity): context removed as well: correction = SMAP - (phi_rem + phi_ctx + base_val)
    df_2019['smap_star_rf_ctx'] = df_2019['smap_sm'] - (phi_rem + phi_ctx + base_val)
    df_2019['d_res_rf_ctx'] = df_2019['smap_star_rf_ctx'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_rf_ctx'] = (df_2019['smap_star_rf_ctx'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_rf_ctx'] = (
        (df_2019['under_wlt_rf_ctx'] == 0) &
        (df_2019['smap_star_rf_ctx'] >= 0.01) &
        (df_2019['smap_star_rf_ctx'] <= 1.0) &
        (~df_2019['smap_star_rf_ctx'].isna())
    )
    
    # Common sample across ALL 5 references simultaneously
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b0'] &
        df_2019['assimilable_b1'] &
        df_2019['assimilable_b2'] &
        df_2019['assimilable_rf'] &
        df_2019['assimilable_rf_ctx']
    )
    
    return df_2019, n_total_2019, n_base_qc


def compute_metrics(df_2019, n_total_2019, n_base_qc):
    """Compute comprehensive metrics on specific and common samples for all methods."""
    print("--- Computing Comparative Stratified Metrics (Specific & Common) ---")
    methods = [
        ('b0', 'd_res_b0', 'assimilable_b0'),
        ('b1', 'd_res_b1', 'assimilable_b1'),
        ('b2', 'd_res_b2', 'assimilable_b2'),
        ('rf', 'd_res_rf', 'assimilable_rf'),
        ('rf_ctx', 'd_res_rf_ctx', 'assimilable_rf_ctx')
    ]
    
    top_textures = ["Loam", "Clay Loam", "Sandy Clay Loam", "Sandy Loam", "Clay"]
    seasons = ["DJF", "MAM", "JJA", "SON"]
    gmia_classes = ["0%", "0-10%", "10-30%", ">30%"]
    
    results = []
    
    def eval_group(sub_df, stratum_type, stratum_name):
        n_sub_base = len(sub_df)
        rec = {
            'stratum_type': stratum_type,
            'stratum_name': stratum_name,
            'total_2019_rows': n_total_2019,
            'base_qc_rows': n_base_qc,
            'stratum_base_qc_rows': n_sub_base,
        }
        
        # Common sample size
        n_com = int(sub_df['assimilable_common'].sum())
        rec['n_common'] = n_com
        rec['pct_common'] = (n_com / n_sub_base * 100.0) if n_sub_base > 0 else 0.0
        
        for m_name, d_col, mask_col in methods:
            sub_spec = sub_df[sub_df[mask_col]]
            n_spec = len(sub_spec)
            pct_spec = (n_spec / n_sub_base * 100.0) if n_sub_base > 0 else 0.0
            bias_spec = float(sub_spec[d_col].mean()) if n_spec > 0 else np.nan
            rmse_spec = float(np.sqrt((sub_spec[d_col] ** 2).mean())) if n_spec > 0 else np.nan
            
            sub_com = sub_df[sub_df['assimilable_common']]
            bias_com = float(sub_com[d_col].mean()) if n_com > 0 else np.nan
            rmse_com = float(np.sqrt((sub_com[d_col] ** 2).mean())) if n_com > 0 else np.nan
            
            # Wilt fraction
            wlt_col = f'under_wlt_{m_name}'
            n_wlt = int(sub_df[wlt_col].sum())
            pct_wlt = (n_wlt / n_sub_base * 100.0) if n_sub_base > 0 else 0.0
            
            rec[f'n_{m_name}_spec'] = n_spec
            rec[f'pct_{m_name}_spec'] = pct_spec
            rec[f'bias_{m_name}_spec'] = bias_spec
            rec[f'rmse_{m_name}_spec'] = rmse_spec
            rec[f'bias_{m_name}_com'] = bias_com
            rec[f'rmse_{m_name}_com'] = rmse_com
            rec[f'n_{m_name}_wlt'] = n_wlt
            rec[f'pct_{m_name}_wlt'] = pct_wlt
            
        results.append(rec)
        
    # 1. Overall
    eval_group(df_2019, "Overall", "Overall")
    # 2. Season
    for s in seasons:
        eval_group(df_2019[df_2019['season'] == s], "Season", s)
    # 3. Texture
    for tex in top_textures:
        eval_group(df_2019[df_2019['texture_name'] == tex], "Texture", tex)
    # 4. GMIA
    for gmia in gmia_classes:
        eval_group(df_2019[df_2019['gmia_class'] == gmia], "GMIA", gmia)
        
    return pd.DataFrame(results)


def compute_anomaly_correlations(df_2019):
    """Compute anomaly correlation with raw SMAP on specific and common assimilable samples."""
    print("--- Computing Anomaly Correlations with Raw SMAP ---")
    methods = [
        ('b0', 'smap_star_b0', 'assimilable_b0'),
        ('b1', 'smap_star_b1', 'assimilable_b1'),
        ('b2', 'smap_star_b2', 'assimilable_b2'),
        ('rf', 'smap_star_rf', 'assimilable_rf'),
        ('rf_ctx', 'smap_star_rf_ctx', 'assimilable_rf_ctx')
    ]
    
    results = {}
    for sample_name, mask in [('specific', None), ('common', df_2019['assimilable_common'])]:
        res_sample = {}
        for m_name, star_col, mask_col in methods:
            sub = df_2019[df_2019[mask_col]] if mask is None else df_2019[mask]
            anom_raw = (sub['smap_sm'] - sub.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')).values
            anom_star = (sub[star_col] - sub.groupby(['pixel_id', 'month'])[star_col].transform('mean')).values
            valid_m = ~(np.isnan(anom_raw) | np.isnan(anom_star))
            if valid_m.sum() > 10:
                r = float(np.corrcoef(anom_raw[valid_m], anom_star[valid_m])[0, 1])
            else:
                r = np.nan
            res_sample[f'r_anom_{m_name}'] = r
            res_sample[f'n_{m_name}'] = len(sub)
        results[sample_name] = res_sample
    return results


def compute_gmia_contrast_with_bootstrap(df_2019):
    """Compute GMIA >30% vs 0% contrast for residual d with 1° block bootstrap 95% CIs and distinct pixel + block counts."""
    print("--- Computing GMIA Contrast (>30% vs 0%) at Fixed Texture and Season ---")
    np.random.seed(42)
    B = 500
    rows = []
    
    target_strata = [
        ("Clay Loam", "All"),
        ("Clay Loam", "DJF"),
        ("Clay Loam", "MAM"),
        ("Clay Loam", "JJA"),
        ("Clay Loam", "SON"),
        ("Loam", "All"),
        ("Loam", "DJF"),
        ("Loam", "MAM"),
        ("Loam", "JJA"),
        ("Loam", "SON"),
    ]
    
    methods = [
        ('b0', 'd_res_b0', 'assimilable_b0'),
        ('b1', 'd_res_b1', 'assimilable_b1'),
        ('b2', 'd_res_b2', 'assimilable_b2'),
        ('rf', 'd_res_rf', 'assimilable_rf'),
        ('rf_ctx', 'd_res_rf_ctx', 'assimilable_rf_ctx')
    ]
    
    for tex, s in target_strata:
        base_sub = df_2019[df_2019['texture_name'] == tex]
        if s != "All":
            base_sub = base_sub[base_sub['season'] == s]
            
        for sample_type in ['specific', 'common']:
            entry = {
                'texture': tex,
                'season': s,
                'sample_type': sample_type
            }
            
            for ref, col_res, mask_col_base in methods:
                mask_col = mask_col_base if sample_type == 'specific' else 'assimilable_common'
                sub_ref = base_sub[base_sub[mask_col] == True]
                sub_0 = sub_ref[sub_ref['gmia_class'] == '0%']
                sub_gt30 = sub_ref[sub_ref['gmia_class'] == '>30%']
                
                n_obs_0 = len(sub_0)
                n_px_0 = int(sub_0['pixel_id'].nunique())
                n_blocks_0 = int(sub_0['spatial_block_id'].nunique())
                n_obs_gt30 = len(sub_gt30)
                n_px_gt30 = int(sub_gt30['pixel_id'].nunique())
                n_blocks_gt30 = int(sub_gt30['spatial_block_id'].nunique())
                
                mean_0 = float(sub_0[col_res].mean()) if n_obs_0 > 0 else np.nan
                mean_gt30 = float(sub_gt30[col_res].mean()) if n_obs_gt30 > 0 else np.nan
                contrast = mean_gt30 - mean_0 if (not np.isnan(mean_0) and not np.isnan(mean_gt30)) else np.nan
                
                sub_both = sub_ref[sub_ref['gmia_class'].isin(['0%', '>30%'])]
                n_blocks_total = int(sub_both['spatial_block_id'].nunique())
                
                if n_obs_0 >= 10 and n_obs_gt30 >= 10 and n_blocks_total > 1:
                    stats = sub_both.groupby(['spatial_block_id', 'gmia_class'])[col_res].agg(['sum', 'count']).unstack('gmia_class').fillna(0)
                    sum_0_arr = stats[('sum', '0%')].values if ('sum', '0%') in stats else np.zeros(len(stats))
                    cnt_0_arr = stats[('count', '0%')].values if ('count', '0%') in stats else np.zeros(len(stats))
                    sum_gt30_arr = stats[('sum', '>30%')].values if ('sum', '>30%') in stats else np.zeros(len(stats))
                    cnt_gt30_arr = stats[('count', '>30%')].values if ('count', '>30%') in stats else np.zeros(len(stats))
                    n_b = len(stats)
                    
                    idx = np.random.choice(n_b, size=(B, n_b), replace=True)
                    c0 = cnt_0_arr[idx].sum(axis=1)
                    cgt = cnt_gt30_arr[idx].sum(axis=1)
                    valid_b = (c0 > 0) & (cgt > 0)
                    if valid_b.sum() > 50:
                        m0 = sum_0_arr[idx].sum(axis=1)[valid_b] / c0[valid_b]
                        mgt = sum_gt30_arr[idx].sum(axis=1)[valid_b] / cgt[valid_b]
                        boot_diffs = mgt - m0
                        ci_low = float(np.percentile(boot_diffs, 2.5))
                        ci_high = float(np.percentile(boot_diffs, 97.5))
                    else:
                        ci_low, ci_high = np.nan, np.nan
                else:
                    ci_low, ci_high = np.nan, np.nan
                    
                entry[f'{ref}_n_obs_0'] = n_obs_0
                entry[f'{ref}_n_px_0'] = n_px_0
                entry[f'{ref}_n_blocks_0'] = n_blocks_0
                entry[f'{ref}_n_obs_gt30'] = n_obs_gt30
                entry[f'{ref}_n_px_gt30'] = n_px_gt30
                entry[f'{ref}_n_blocks_gt30'] = n_blocks_gt30
                entry[f'{ref}_n_blocks_total'] = n_blocks_total
                entry[f'{ref}_contrast'] = contrast
                entry[f'{ref}_ci_low'] = ci_low
                entry[f'{ref}_ci_high'] = ci_high
                
            rows.append(entry)
            
    return pd.DataFrame(rows)


def generate_rf_report(df_metrics, df_cv, df_shap_stab, df_contrast, anom_corr, n_total_2019, n_base_qc, base_val, n_train_rows, n_train_pixels):
    """Generate Phase 2b Markdown report with zero hardcoded numbers and zero interpretive words."""
    print("--- Generating Phase 2b Report (Zero Interpretive Phrasing) ---")
    out_md = os.path.join(REPORTS_DIR, "phase2b_rf_report.md")
    
    glob = df_metrics[df_metrics['stratum_type'] == 'Overall'].iloc[0]
    
    # CV table rows
    cv_rows = []
    for _, r in df_cv.iterrows():
        cv_rows.append(
            f"| {int(r['spatial_fold'])} | {int(r['n_train']):,} | {int(r['n_val']):,} | "
            f"{r['oof_rmse']:.4f} | {r['oof_mae']:.4f} | {r['oof_bias']:+.4f} | {r['oof_r2']:.4f} |"
        )
    cv_table = chr(10).join(cv_rows)
    
    # SHAP stability table rows
    shap_rows = []
    for _, r in df_shap_stab.iterrows():
        shap_rows.append(
            f"| Plis {int(r['spatial_fold'])} | {r['mean_abs_removable']:.4f} ({r['share_removable_pct']:.1f}%) | "
            f"{r['mean_abs_retained']:.4f} ({r['share_retained_pct']:.1f}%) | "
            f"{r['mean_abs_context']:.4f} ({r['share_context_pct']:.1f}%) |"
        )
    shap_table = chr(10).join(shap_rows)
    
    # Mean attribution row
    mean_rem = df_shap_stab['share_removable_pct'].mean()
    mean_ret = df_shap_stab['share_retained_pct'].mean()
    mean_ctx = df_shap_stab['share_context_pct'].mean()
    shap_mean_row = f"| **Moyenne 5 plis** | **{mean_rem:.1f}%** | **{mean_ret:.1f}%** | **{mean_ctx:.1f}%** |"
    
    # Population table rows
    pop_rows = [
        f"| Spécifique | **B0 (NoCDF)** | {glob['n_b0_spec']:,} | {glob['pct_b0_spec']:.2f}% | {glob['n_b0_spec']/n_total_2019*100:.2f}% |",
        f"| Spécifique | **B1 (Climatologie)** | {glob['n_b1_spec']:,} | {glob['pct_b1_spec']:.2f}% | {glob['n_b1_spec']/n_total_2019*100:.2f}% |",
        f"| Spécifique | **B2 (CDF v2)** | {glob['n_b2_spec']:,} | {glob['pct_b2_spec']:.2f}% | {glob['n_b2_spec']/n_total_2019*100:.2f}% |",
        f"| Spécifique | **RF (SAR main)** | {glob['n_rf_spec']:,} | {glob['pct_rf_spec']:.2f}% | {glob['n_rf_spec']/n_total_2019*100:.2f}% |",
        f"| Spécifique | **RF-ctx (Sensibilité)** | {glob['n_rf_ctx_spec']:,} | {glob['pct_rf_ctx_spec']:.2f}% | {glob['n_rf_ctx_spec']/n_total_2019*100:.2f}% |",
        f"| Commun | **Commun (5 méthodes)** | {glob['n_common']:,} | {glob['pct_common']:.2f}% | {glob['n_common']/n_total_2019*100:.2f}% |"
    ]
    pop_table = chr(10).join(pop_rows)
    
    # Global metrics table
    glob_table_rows = [
        f"| **B0 (NoCDF)** | {glob['bias_b0_spec']:+.4f} m³/m³ | {glob['rmse_b0_spec']:.4f} m³/m³ | {glob['bias_b0_com']:+.4f} m³/m³ | {glob['rmse_b0_com']:.4f} m³/m³ |",
        f"| **B1 (Climatologie)** | {glob['bias_b1_spec']:+.4f} m³/m³ | {glob['rmse_b1_spec']:.4f} m³/m³ | {glob['bias_b1_com']:+.4f} m³/m³ | {glob['rmse_b1_com']:.4f} m³/m³ |",
        f"| **B2 (CDF v2)** | {glob['bias_b2_spec']:+.4f} m³/m³ | {glob['rmse_b2_spec']:.4f} m³/m³ | {glob['bias_b2_com']:+.4f} m³/m³ | {glob['rmse_b2_com']:.4f} m³/m³ |",
        f"| **RF (SAR main)** | {glob['bias_rf_spec']:+.4f} m³/m³ | {glob['rmse_rf_spec']:.4f} m³/m³ | {glob['bias_rf_com']:+.4f} m³/m³ | {glob['rmse_rf_com']:.4f} m³/m³ |",
        f"| **RF-ctx (Sensibilité)** | {glob['bias_rf_ctx_spec']:+.4f} m³/m³ | {glob['rmse_rf_ctx_spec']:.4f} m³/m³ | {glob['bias_rf_ctx_com']:+.4f} m³/m³ | {glob['rmse_rf_ctx_com']:.4f} m³/m³ |"
    ]
    glob_metrics_table = chr(10).join(glob_table_rows)
    
    # Seasonal RMSE table
    season_rows = []
    for s in ["DJF", "MAM", "JJA", "SON"]:
        s_row = df_metrics[df_metrics['stratum_name'] == s].iloc[0]
        season_rows.append(
            f"| {s} | {s_row['rmse_b0_com']:.4f} | {s_row['rmse_b1_com']:.4f} | {s_row['rmse_b2_com']:.4f} | "
            f"{s_row['rmse_rf_com']:.4f} | {s_row['rmse_rf_ctx_com']:.4f} | {s_row['n_common']:,} |"
        )
    season_table = chr(10).join(season_rows)
    
    # Anomaly correlation table
    an_sp = anom_corr['specific']
    an_co = anom_corr['common']
    anom_rows = [
        f"| **B0 (NoCDF)** | {an_sp['r_anom_b0']:.4f} | {an_sp['n_b0']:,} | {an_co['r_anom_b0']:.4f} | {an_co['n_b0']:,} |",
        f"| **B1 (Climatologie)** | {an_sp['r_anom_b1']:.4f} | {an_sp['n_b1']:,} | {an_co['r_anom_b1']:.4f} | {an_co['n_b1']:,} |",
        f"| **B2 (CDF v2)** | {an_sp['r_anom_b2']:.4f} | {an_sp['n_b2']:,} | {an_co['r_anom_b2']:.4f} | {an_co['n_b2']:,} |",
        f"| **RF (SAR main)** | {an_sp['r_anom_rf']:.4f} | {an_sp['n_rf']:,} | {an_co['r_anom_rf']:.4f} | {an_co['n_rf']:,} |",
        f"| **RF-ctx (Sensibilité)** | {an_sp['r_anom_rf_ctx']:.4f} | {an_sp['n_rf_ctx']:,} | {an_co['r_anom_rf_ctx']:.4f} | {an_co['n_rf_ctx']:,} |"
    ]
    anom_table = chr(10).join(anom_rows)
    
    # Contrast table (Clay Loam across seasons)
    contrast_rows = []
    for _, c_row in df_contrast[df_contrast['texture'] == 'Clay Loam'].iterrows():
        def fmt_ci(c, low, high):
            if np.isnan(c) or np.isnan(low) or np.isnan(high):
                return "N/A"
            return f"{c:+.4f} [{low:+.4f}, {high:+.4f}]"

        contrast_rows.append(
            f"| {c_row['season']} | {c_row['sample_type']} | "
            f"{fmt_ci(c_row['b0_contrast'], c_row['b0_ci_low'], c_row['b0_ci_high'])} | "
            f"{fmt_ci(c_row['b1_contrast'], c_row['b1_ci_low'], c_row['b1_ci_high'])} | "
            f"{fmt_ci(c_row['b2_contrast'], c_row['b2_ci_low'], c_row['b2_ci_high'])} | "
            f"{fmt_ci(c_row['rf_contrast'], c_row['rf_ci_low'], c_row['rf_ci_high'])} | "
            f"{c_row['rf_n_px_gt30']}/{c_row['rf_n_px_0']} px, {c_row['rf_n_blocks_gt30']}/{c_row['rf_n_blocks_0']} blk |"
        )
    contrast_table = chr(10).join(contrast_rows)
    
    # Seasonal under wilt table
    wlt_rows = []
    for s in ["DJF", "MAM", "JJA", "SON"]:
        s_row = df_metrics[df_metrics['stratum_name'] == s].iloc[0]
        wlt_rows.append(
            f"| {s} | {s_row['pct_b0_wlt']:.2f}% | {s_row['pct_b1_wlt']:.2f}% | "
            f"{s_row['pct_b2_wlt']:.2f}% | {s_row['pct_rf_wlt']:.2f}% | {s_row['pct_rf_ctx_wlt']:.2f}% |"
        )
    wlt_table = chr(10).join(wlt_rows)
    
    latex_contrast = r"$\Delta d_{\text{res}} = \bar{d}_{\text{res}}(>30\%) - \bar{d}_{\text{res}}(0\%)$"
    latex_wlt = r"$\text{SM}^* < \theta_{\text{wlt}} + 0.02$"
    latex_anom = r"$r_{\text{anom}}$"
    
    report_md = f"""# Rapport d'Évaluation Phase 2b : Premier Modèle d'Apprentissage (Random Forest)

**Date d'exécution** : {time.strftime('%Y-%m-%d %H:%M:%S')}  
**Période d'entraînement** : 2016–2018 (`split == 'train'`).  
**Période d'évaluation** : 2019 (`split == 'val'`).  
**Statut de l'année 2020** : 2020 reste strictement scellée.  
**Sous-échantillonnage d'entraînement** : Échantillonnage stratifié uniforme de 20 observations par pixel sur les observations conformes aux règles QC hors flétrissement. Effectif d'entraînement = {n_train_rows:,} observations réparties uniformément sur {n_train_pixels:,} mailles continentales.  
**Valeur de base (phi_0)** : {base_val:+.6f} m³/m³.

---

## 1. Validation Croisée Spatiale (5 Plis Équilibrés de 1°)

| Pli de validation | Effectif train | Effectif val | RMSE hors-pli | MAE hors-pli | Biais hors-pli | R² hors-pli |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{cv_table}

---

## 2. Attributions SHAP par Groupe et Stabilité Spatiale (Jalon H4)

Attributions SHAP moyennes calculées sur les échantillons de validation de chaque pli :

| Pli spatial | Groupe Removable (à retirer) | Groupe Retained (à conserver) | Groupe Context (analyse sensibilité) |
|:---|:---:|:---:|:---:|
{shap_table}
{shap_mean_row}

---

## 3. Population d'Évaluation 2019 et Observations Assimilables Après Recalcul du Flétrissement

Socle d'observations 2019 passant les règles QC hors flétrissement : {n_base_qc:,} observations ({n_base_qc/n_total_2019*100:.2f}% des {n_total_2019:,} lignes brutes).

| Échantillon | Méthode | Observations assimilables | Part du socle hors flétrissement ({n_base_qc:,} obs) | Part du total 2019 ({n_total_2019:,} obs) |
|:---|:---|:---:|:---:|:---:|
{pop_table}

---

## 4. Métriques Globales Comparatives : Biais et RMSE Résiduels

| Méthode | Biais résiduel (Échantillon spécifique) | RMSE résiduelle (Échantillon spécifique) | Biais résiduel (Échantillon commun) | RMSE résiduelle (Échantillon commun) |
|:---|:---:|:---:|:---:|:---:|
{glob_metrics_table}

---

## 5. Performance Saisonnière sur l'Échantillon Commun (RMSE Résiduelle)

| Saison | RMSE B0 (NoCDF) | RMSE B1 (Climatologie) | RMSE B2 (CDF v2) | RMSE RF (SAR main) | RMSE RF-ctx (Sensibilité) | Observations communes |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
{season_table}

---

## 6. Corrélation d'Anomalies avec SMAP Brut ({latex_anom})

| Méthode | {latex_anom} spécifique | N spécifique | {latex_anom} commun | N commun |
|:---|:---:|:---:|:---:|:---:|
{anom_table}

---

## 7. Préservation du Signal d'Irrigation : Contraste GMIA >30 % vs 0 % en Clay Loam

Contraste {latex_contrast} avec intervalle de confiance à 95% par bootstrap par blocs de 1° (500 itérations) :

| Saison | Échantillon | Contraste B0 | Contraste B1 | Contraste B2 | Contraste RF (SAR) | Effectif (>30% / 0%) |
|:---|:---|:---:|:---:|:---:|:---:|:---:|
{contrast_table}

---

## 8. Part des Valeurs Sous le Seuil de Flétrissement par Saison ({latex_wlt})

| Saison | B0 (NoCDF) | B1 (Climatologie) | B2 (CDF v2) | RF (SAR main) | RF-ctx (Sensibilité) |
|:---|:---:|:---:|:---:|:---:|
{wlt_table}
## 9. Fichiers et Tables Produites

- **Table de métriques détaillées** : `reports/phase2b_rf_metrics.csv`
- **Table de stabilité SHAP 5 plis** : `reports/phase2b_rf_shap_stability.csv`
- **Table de contrastes GMIA** : `reports/phase2b_rf_gmia_contrast.csv`
"""
    with open(out_md, 'w') as f:
        f.write(report_md)
    print(f"Generated report: {out_md}")


def main():
    print("=== STARTING PHASE 2b RANDOM FOREST MODELING ===")
    t_start = time.time()
    
    # 1. Feature groups
    cols_rem, cols_ret, cols_ctx = get_feature_groups()
    feature_cols = cols_rem + cols_ret + cols_ctx
    print(f"Predictor count: {len(feature_cols)} (Removable: {len(cols_rem)}, Retained: {len(cols_ret)}, Context: {len(cols_ctx)})")
    
    # 2. Static parameters
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]
        
    # 3. Training data
    train_df = load_training_data(feature_cols)
    
    # 4. 5-Fold spatial CV and SHAP stability
    df_cv, df_shap_stab = run_spatial_cross_validation(train_df, feature_cols, cols_rem, cols_ret, cols_ctx)
    csv_shap_stab = os.path.join(REPORTS_DIR, "phase2b_rf_shap_stability.csv")
    df_shap_stab.to_csv(csv_shap_stab, index=False)
    print(f"Saved SHAP stability CSV: {csv_shap_stab}")
    csv_cv = os.path.join(REPORTS_DIR, "phase2b_rf_cv_results.csv")
    df_cv.to_csv(csv_cv, index=False)
    print(f"Saved CV results CSV: {csv_cv}")
    
    # 5. Full Model training
    rf_full, explainer, base_val = train_full_model(train_df, feature_cols)
    
    # 6. Baseline assets
    lut, climatology, pixel_fallback, global_fallback = load_lut_and_climatology(landmask)
    
    # 7. Evaluate on 2019
    df_2019, n_total_2019, n_base_qc = evaluate_phase2b_on_2019(
        rf_full, explainer, base_val, feature_cols, cols_rem, cols_ret, cols_ctx,
        lut, climatology, pixel_fallback, global_fallback, landmask, soil_wlt, dom_soil
    )
    
    # 8. Metrics calculation
    df_metrics = compute_metrics(df_2019, n_total_2019, n_base_qc)
    csv_metrics = os.path.join(REPORTS_DIR, "phase2b_rf_metrics.csv")
    df_metrics.to_csv(csv_metrics, index=False)
    print(f"Saved metrics CSV: {csv_metrics}")
    
    # 9. Anomaly correlations
    anom_corr = compute_anomaly_correlations(df_2019)
    
    # 10. GMIA contrast with bootstrap
    df_contrast = compute_gmia_contrast_with_bootstrap(df_2019)
    csv_contrast = os.path.join(REPORTS_DIR, "phase2b_rf_gmia_contrast.csv")
    df_contrast.to_csv(csv_contrast, index=False)
    print(f"Saved contrast CSV: {csv_contrast}")
    
    # 11. Generate Markdown Report
    generate_rf_report(df_metrics, df_cv, df_shap_stab, df_contrast, anom_corr, n_total_2019, n_base_qc, base_val, len(train_df), int(train_df['pixel_id'].nunique()))
    
    print(f"=== PHASE 2b COMPLETED IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
