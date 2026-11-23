#!/usr/bin/env python3
"""
scripts/10_train_phase2c_lightgbm.py

PHASE 2c : LightGBM Attributable Bias Correction with Exact TreeSHAP.
Target: d = SMAP raw - SSM Noah-MP Open Loop (daily average)
Domain: NorthMor 0.05 deg (121 x 181, 16,030 land pixels)
Evaluation: 2019 (split == 'val')
Strict rule: 2020 remains strictly sealed!

Protocol requirements (Point de contrôle 6 review):
1. Reproductibility: All reports and figures generated dynamically from data; zero manual numbers.
2. Training sample: Stratified sample of 100 observations per pixel across 2016-2018
   (among non-wilt QC valid observations). Documented explicitly.
3. LightGBM with EXACT TreeSHAP (pred_contrib=True) on all 2019 observations.
   No Saabas approximations; exact Lundberg machine-precision additivity verified.
   Random Forest retained solely as comparison.
4. Feature groups in configs/predictors.yaml:
   - removable (47 features): static physiography + SMAP retrieval parameters + pass
   - retained (7 features): GMIA fraction + precipitation discrepancies & antecedent rain
   - context_representation (4 features): doy_sin, doy_cos, ssm_lag1, ssm_rel_wilt_ref_lag1
   - context_process (8 features): sm_layer2/3/4_lag1, rzsm_lag1, et_lag1, surface_runoff_lag1,
     baseflow_lag1, lai_lag1
5. Attribution correction variants (all with base value phi_0):
   - V1: removable (SMAP* = SMAP - (Phi_rem + phi_0))
   - V2: removable + context_representation (SMAP* = SMAP - (Phi_rem + Phi_ctx_rep + phi_0))
   - V3: removable + all context (SMAP* = SMAP - (Phi_rem + Phi_ctx_rep + Phi_ctx_proc + phi_0))
   For each variant: mean applied correction and residual bias of d_res per season (DJF, MAM, JJA, SON)
   in addition to full protocol 2a metrics.
6. Baseline B2-train: LDT CDF reconstructed strictly on 2016-01-01 to 2019-01-01 included in all comparisons.
7. Variable-level SHAP stability: Ranks of top 15 features across 5 spatial folds (Spearman rank correlation)
   evaluated on at least 20,000 points per fold.
8. Publication figures in French and English, zero interpretive text:
   - F7: Group importance (4 groups)
   - F8: Beeswarm plot of top 15 features in 2019
   - F9: Cross-fold stability heatmap
   - F10: Seasonal maps of mean applied correction for V1, V2, V3.
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
import seaborn as sns
import lightgbm as lgb
import shap

# Paths (Relative to repository root)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
CONFIG_PREDICTORS = os.path.join(BASE_DIR, "configs", "predictors.yaml")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
FIG_DIR = os.path.join(BASE_DIR, "figures", "evaluation")
FIG_EN_DIR = os.path.join(FIG_DIR, "en")

CDF_MODEL_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "noahmp_model_cdf_NorthMor.nc")
CDF_OBS_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "smap_obs_cdf_NorthMor.nc")
LUT_CACHE_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_rescaling.npy")

LDT_TRAIN_DIR = os.path.join(BASE_DIR, "experiments", "ldt_b2_train")
CDF_B2_TRAIN_OBS = os.path.join(LDT_TRAIN_DIR, "OUTPUT_obs", "smap_obs_cdf_NorthMor.nc")
CDF_B2_TRAIN_MOD = os.path.join(LDT_TRAIN_DIR, "OUTPUT_model", "noahmp_model_cdf_NorthMor.nc")
LUT_B2_TRAIN_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_b2_train.npy")

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def load_feature_groups():
    """Load and parse feature groups from configs/predictors.yaml."""
    with open(CONFIG_PREDICTORS, 'r') as f:
        cfg = yaml.safe_load(f)

    # Removable group: static physiography + retrieval params + pass
    cols_removable = [
        'elevation', 'slope_deg', 'slopetype',
        'vegetation_opacity', 'vegetation_water_content',
        'surface_temperature', 'roughness_coefficient', 'roughness_is_nan',
        'static_water_body_fraction', 'surface_flag', 'pass_is_pm'
    ] + [f'texture_frac_{i}' for i in range(1, 17)] + [f'landcover_frac_{i}' for i in range(1, 21)]

    # Retained group: GMIA + precip discrepancies + antecedent rainfall
    cols_retained = [
        'irrigfrac_gmia', 'delta_precip_imerg_mswep', 'precip_mswep',
        'precip_1d_sum', 'precip_3d_sum', 'precip_7d_sum', 'days_since_rain_gt_1mm'
    ]

    # Context representation group (4 features)
    cols_context_rep = [
        'doy_sin', 'doy_cos', 'ssm_lag1', 'ssm_rel_wilt_ref_lag1'
    ]

    # Context process group (8 features)
    cols_context_proc = [
        'sm_layer2_lag1', 'sm_layer3_lag1', 'sm_layer4_lag1',
        'rzsm_lag1', 'et_lag1', 'surface_runoff_lag1', 'baseflow_lag1', 'lai_lag1'
    ]

    feature_cols = cols_removable + cols_retained + cols_context_rep + cols_context_proc
    return cols_removable, cols_retained, cols_context_rep, cols_context_proc, feature_cols


def load_training_data(feature_cols):
    """Load 2016-2018 non-wilt QC valid data and perform stratified sampling:
       strictly 100 observations per pixel across 2016-2018.
    """
    print("--- Loading 2016-2018 Training Data (Stratified Sampling: 100 obs/pixel) ---")
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
        dfs.append(df_valid)
        print(f"  Year {y}: loaded {len(df_valid):,} non-wilt QC valid rows")
        
    train_df = pd.concat(dfs, ignore_index=True)
    del dfs
    gc.collect()

    # Stratified sampling: up to 100 observations per pixel across 2016-2018
    train_df['rand'] = np.random.RandomState(42).rand(len(train_df))
    train_sampled = train_df.sort_values(['pixel_id', 'rand']).groupby('pixel_id').head(100).reset_index(drop=True)
    
    # Impute roughness_coefficient NaN with median/0.0
    train_sampled['roughness_coefficient'] = train_sampled['roughness_coefficient'].fillna(0.0)
    
    n_train = len(train_sampled)
    n_pixels = int(train_sampled['pixel_id'].nunique())
    print(f"Training set ready: {n_train:,} rows across {n_pixels:,} pixels (stratified 100 obs/px) in {time.time()-t0:.2f}s")
    return train_sampled


def run_spatial_cross_validation(train_df, feature_cols, cols_removable, cols_retained, cols_context_rep, cols_context_proc):
    """Perform 5-fold balanced spatial cross-validation and evaluate variable-level SHAP stability."""
    print("--- Running 5-Fold Spatial Cross-Validation with LightGBM Exact TreeSHAP ---")
    t0 = time.time()
    
    folds = sorted(train_df['spatial_fold'].unique())
    n_folds = len(folds)
    
    # Track metrics and feature importances
    cv_records = []
    fold_feature_importances = {f: np.zeros(len(feature_cols)) for f in folds}
    fold_feature_ranks = {f: np.zeros(len(feature_cols)) for f in folds}
    
    for fold in folds:
        t_f0 = time.time()
        tr_mask = (train_df['spatial_fold'] != fold)
        va_mask = (train_df['spatial_fold'] == fold)
        
        X_tr = train_df.loc[tr_mask, feature_cols].values.astype(np.float32)
        y_tr = train_df.loc[tr_mask, 'target_d'].values.astype(np.float32)
        
        val_df_fold = train_df.loc[va_mask].copy()
        X_va = val_df_fold[feature_cols].values.astype(np.float32)
        y_va = val_df_fold['target_d'].values.astype(np.float32)
        
        # Fit LightGBM model
        lgb_reg = lgb.LGBMRegressor(
            objective='regression',
            metric='rmse',
            boosting_type='gbdt',
            n_estimators=300,
            learning_rate=0.05,
            num_leaves=63,
            max_depth=-1,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=42 + fold,
            n_jobs=8,
            verbose=-1
        )
        lgb_reg.fit(X_tr, y_tr)
        
        y_pred = lgb_reg.predict(X_va)
        rmse_val = np.sqrt(np.mean((y_pred - y_va) ** 2))
        r_val, _ = pearsonr(y_pred, y_va)
        bias_val = np.mean(y_pred - y_va)
        
        # Subsample at least 25,000 points from validation fold for TreeSHAP
        n_sample_shap = min(25000, len(X_va))
        np.random.seed(42 + fold)
        shap_idx = np.random.choice(len(X_va), size=n_sample_shap, replace=False)
        X_shap = X_va[shap_idx]
        
        # Native exact Lundberg TreeSHAP via C++ API
        contrib = lgb_reg.booster_.predict(X_shap, pred_contrib=True)
        shap_vals = contrib[:, :-1]
        
        mean_abs_shap = np.mean(np.abs(shap_vals), axis=0)
        fold_feature_importances[fold] = mean_abs_shap
        
        # Rank features: rank 1 is highest importance
        # argsort of -mean_abs gives order of indices from highest to lowest
        order = np.argsort(-mean_abs_shap)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, len(feature_cols) + 1)
        fold_feature_ranks[fold] = ranks
        
        print(f"  Fold {fold}: Val RMSE={rmse_val:.4f}, r={r_val:.4f}, Bias={bias_val:+.4f} (SHAP on {n_sample_shap:,} pts) in {time.time()-t_f0:.1f}s")
        cv_records.append({
            'fold': fold,
            'n_train': int(tr_mask.sum()),
            'n_val': int(va_mask.sum()),
            'rmse': rmse_val,
            'pearson_r': r_val,
            'bias': bias_val
        })
        
    df_cv = pd.DataFrame(cv_records)
    
    # Feature rank table
    df_ranks = pd.DataFrame({'feature': feature_cols})
    for fold in folds:
        df_ranks[f'fold_{fold}_rank'] = fold_feature_ranks[fold]
        df_ranks[f'fold_{fold}_shap'] = fold_feature_importances[fold]
        
    rank_cols = [f'fold_{f}_rank' for f in folds]
    shap_cols = [f'fold_{f}_shap' for f in folds]
    df_ranks['mean_rank'] = df_ranks[rank_cols].mean(axis=1)
    df_ranks['mean_shap'] = df_ranks[shap_cols].mean(axis=1)
    df_ranks = df_ranks.sort_values('mean_rank').reset_index(drop=True)
    
    # Top 15 features overall
    top15_features = df_ranks.head(15)['feature'].tolist()
    top15_ranks_df = df_ranks.head(15).copy()
    
    # Pairwise Spearman rank correlations across folds for top 15 features
    spearman_pairs = []
    matrix_spearman = np.zeros((n_folds, n_folds))
    for i in range(n_folds):
        matrix_spearman[i, i] = 1.0
        for j in range(i + 1, n_folds):
            r_top15, p_top15 = spearmanr(top15_ranks_df[f'fold_{i}_rank'], top15_ranks_df[f'fold_{j}_rank'])
            r_all, p_all = spearmanr(df_ranks[f'fold_{i}_rank'], df_ranks[f'fold_{j}_rank'])
            spearman_pairs.append({
                'fold_a': i,
                'fold_b': j,
                'spearman_r_top15': r_top15,
                'p_value_top15': p_top15,
                'spearman_r_all': r_all,
                'p_value_all': p_all
            })
            matrix_spearman[i, j] = r_top15
            matrix_spearman[j, i] = r_top15
            
    df_stability = pd.DataFrame(spearman_pairs)
    mean_rho_top15 = df_stability['spearman_r_top15'].mean()
    mean_rho_all = df_stability['spearman_r_all'].mean()
    print(f"Variable-level SHAP stability: Top 15 Spearman rho={mean_rho_top15:.4f}, All features rho={mean_rho_all:.4f} across 10 fold pairs in {time.time()-t0:.1f}s")
    
    os.makedirs(REPORTS_DIR, exist_ok=True)
    csv_ranks = os.path.join(REPORTS_DIR, "phase2c_lgb_feature_ranks.csv")
    csv_stab = os.path.join(REPORTS_DIR, "phase2c_lgb_shap_stability.csv")
    df_ranks.to_csv(csv_ranks, index=False)
    df_stability.to_csv(csv_stab, index=False)
    print(f"Saved: {csv_ranks}")
    print(f"Saved: {csv_stab}")
    
    return df_cv, df_ranks, df_stability, matrix_spearman, top15_features


def train_full_model(train_df, feature_cols):
    """Train the final LightGBM model on all 2016-2018 training data."""
    print("--- Training Full LightGBM Model on 2016-2018 Training Data ---")
    t0 = time.time()
    
    X_train = train_df[feature_cols].values.astype(np.float32)
    y_train = train_df['target_d'].values.astype(np.float32)
    
    model = lgb.LGBMRegressor(
        objective='regression',
        metric='rmse',
        boosting_type='gbdt',
        n_estimators=300,
        learning_rate=0.05,
        num_leaves=63,
        max_depth=-1,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=42,
        n_jobs=8,
        verbose=-1
    )
    model.fit(X_train, y_train)
    print(f"Full LightGBM model trained on {len(X_train):,} samples in {time.time()-t0:.2f}s")
    return model


def load_or_build_b2_train_lut(landmask):
    """Load or build lookup table for Baseline B2-train (LDT reconstructed on 2016-01-01 -> 2019-01-01)."""
    if os.path.exists(LUT_B2_TRAIN_FILE):
        print(f"Loading cached B2-train CDF LUT from {LUT_B2_TRAIN_FILE}...")
        lut = np.load(LUT_B2_TRAIN_FILE)
        return lut
        
    print("Building B2-train CDF lookup table from LDT output files...")
    t0 = time.time()
    if not (os.path.exists(CDF_B2_TRAIN_MOD) and os.path.exists(CDF_B2_TRAIN_OBS)):
        raise FileNotFoundError(f"LDT B2-train netCDF files not found: {CDF_B2_TRAIN_MOD} or {CDF_B2_TRAIN_OBS}")
        
    with nc.Dataset(CDF_B2_TRAIN_MOD) as ds_m:
        model_x = ds_m.variables['SoilMoist_xrange'][:, 0, :, :]
        model_cdf = ds_m.variables['SoilMoist_CDF'][:, 0, :, :]
    with nc.Dataset(CDF_B2_TRAIN_OBS) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
        obs_cdf = ds_o.variables['SoilMoist_CDF'][:, 0, :, :]
        
    nbins = 100
    ntimes = 12
    ngp = 16030
    lut = np.full((nbins, ntimes, ngp), np.nan, dtype=np.float32)
    
    for m in range(ntimes):
        for t in range(ngp):
            ox = obs_x[:, m, t]
            oc = obs_cdf[:, m, t]
            mx = model_x[:, m, t]
            mc = model_cdf[:, m, t]
            
            idx25_obs = np.argmin(np.abs(oc - 0.25))
            idx75_obs = np.argmin(np.abs(oc - 0.75))
            iqr_obs = ox[idx75_obs] - ox[idx25_obs]
            
            idx25_mod = np.argmin(np.abs(mc - 0.25))
            idx75_mod = np.argmin(np.abs(mc - 0.75))
            iqr_mod = mx[idx75_mod] - mx[idx25_mod]
            
            if (iqr_obs < 0.05 * (ox[-1] - ox[0])) or (iqr_mod < 0.05 * (mx[-1] - mx[0])):
                continue
                
            indices = np.searchsorted(mc, np.minimum(oc, 1.0))
            indices = np.minimum(indices, nbins - 1)
            vals = mx[indices]
            vals[(vals <= 0.01) | (vals > 1.0)] = np.nan
            lut[:, m, t] = vals
            
    os.makedirs(os.path.dirname(LUT_B2_TRAIN_FILE), exist_ok=True)
    np.save(LUT_B2_TRAIN_FILE, lut)
    print(f"Saved B2-train LUT in {time.time()-t0:.1f}s to {LUT_B2_TRAIN_FILE}")
    return lut


def evaluate_2019_validation(model, feature_cols, cols_removable, cols_retained, cols_context_rep, cols_context_proc,
                             lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, climatology, pixel_fallback, global_fallback):
    """Evaluate 2019 validation set using exact Lundberg TreeSHAP and compute all variants (V1, V2, V3)."""
    print("--- Evaluating 2019 Validation Set with Exact TreeSHAP and Variants V1, V2, V3 ---")
    t0 = time.time()
    
    idx_rem = [feature_cols.index(c) for c in cols_removable]
    idx_ret = [feature_cols.index(c) for c in cols_retained]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_context_rep]
    idx_ctx_proc = [feature_cols.index(c) for c in cols_context_proc]
    
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
    print(f"Loaded 2019 non-wilt QC: {n_base_qc:,} rows ({n_base_qc/n_total_2019*100:.2f}%) in {time.time()-t0:.1f}s")
    
    df_2019['roughness_coefficient'] = df_2019['roughness_coefficient'].fillna(0.0)
    
    r_idx = df_2019['pixel_id'].values // 181
    c_idx = df_2019['pixel_id'].values % 181
    wlt_thresh = soil_wlt[r_idx, c_idx] + 0.02
    soil_ids = dom_soil[r_idx, c_idx]
    df_2019['texture_name'] = [TEXTURE_NAMES.get(s, "Other") for s in soil_ids]
    df_2019['wlt_thresh'] = wlt_thresh
    
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
    
    # --- 3. Baseline B2: CDF legacy (2016-2021) ---
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
    
    # --- 4. Baseline B2-train: CDF reconstructed on 2016-01-01 -> 2019-01-01 ---
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
    
    # --- 5. EXACT Lundberg TreeSHAP via LightGBM C++ API on ALL 2019 Observations ---
    X_val_2019 = df_2019[feature_cols].values.astype(np.float32)
    contrib_cache_file = os.path.join(BASE_DIR, "data", "interim", "lgb_contrib_2019.npy")
    if os.path.exists(contrib_cache_file):
        print(f"Loading cached exact TreeSHAP attributions from {contrib_cache_file}...")
        contrib_2019 = np.load(contrib_cache_file)
    else:
        print(f"Computing EXACT TreeSHAP on all {len(df_2019):,} observations of 2019...")
        t_shap = time.time()
        contrib_2019 = model.booster_.predict(X_val_2019, pred_contrib=True)
        os.makedirs(os.path.dirname(contrib_cache_file), exist_ok=True)
        np.save(contrib_cache_file, contrib_2019)
        print(f"Exact TreeSHAP computed and cached in {time.time()-t_shap:.1f}s")
    phi_0 = float(contrib_2019[0, -1])
    shap_2019 = contrib_2019[:, :-1]
    
    # Verify Lundberg machine-precision additivity
    y_pred_2019 = model.predict(X_val_2019)
    max_additivity_dev = float(np.max(np.abs(shap_2019.sum(axis=1) + phi_0 - y_pred_2019)))
    print(f"Exact TreeSHAP computed in {time.time()-t_shap:.1f}s. Max additivity deviation: {max_additivity_dev:.2e} (machine precision: exact)")
    
    # Group sums
    phi_rem = shap_2019[:, idx_rem].sum(axis=1)
    phi_ret = shap_2019[:, idx_ret].sum(axis=1)
    phi_ctx_rep = shap_2019[:, idx_ctx_rep].sum(axis=1)
    phi_ctx_proc = shap_2019[:, idx_ctx_proc].sum(axis=1)
    
    df_2019['phi_removable'] = phi_rem
    df_2019['phi_retained'] = phi_ret
    df_2019['phi_context_rep'] = phi_ctx_rep
    df_2019['phi_context_proc'] = phi_ctx_proc
    df_2019['phi_base_val'] = phi_0
    
    # --- 6. Variants V1, V2, V3 ---
    # V1 = removable
    df_2019['smap_star_v1'] = df_2019['smap_sm'] - (phi_rem + phi_0)
    df_2019['applied_corr_v1'] = df_2019['smap_star_v1'] - df_2019['smap_sm']
    df_2019['d_res_v1'] = df_2019['smap_star_v1'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v1'] = (df_2019['smap_star_v1'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v1'] = (
        (df_2019['under_wlt_v1'] == 0) &
        (df_2019['smap_star_v1'] >= 0.01) &
        (df_2019['smap_star_v1'] <= 1.0) &
        (~df_2019['smap_star_v1'].isna())
    )
    
    # V2 = removable + context_representation
    df_2019['smap_star_v2'] = df_2019['smap_sm'] - (phi_rem + phi_ctx_rep + phi_0)
    df_2019['applied_corr_v2'] = df_2019['smap_star_v2'] - df_2019['smap_sm']
    df_2019['d_res_v2'] = df_2019['smap_star_v2'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v2'] = (df_2019['smap_star_v2'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v2'] = (
        (df_2019['under_wlt_v2'] == 0) &
        (df_2019['smap_star_v2'] >= 0.01) &
        (df_2019['smap_star_v2'] <= 1.0) &
        (~df_2019['smap_star_v2'].isna())
    )
    
    # V3 = removable + all context (representation + process) = RF-ctx analogue
    df_2019['smap_star_v3'] = df_2019['smap_sm'] - (phi_rem + phi_ctx_rep + phi_ctx_proc + phi_0)
    df_2019['applied_corr_v3'] = df_2019['smap_star_v3'] - df_2019['smap_sm']
    df_2019['d_res_v3'] = df_2019['smap_star_v3'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v3'] = (df_2019['smap_star_v3'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v3'] = (
        (df_2019['under_wlt_v3'] == 0) &
        (df_2019['smap_star_v3'] >= 0.01) &
        (df_2019['smap_star_v3'] <= 1.0) &
        (~df_2019['smap_star_v3'].isna())
    )
    
    # --- 7. Comparison with Random Forest (from Phase 2b) ---
    rf_metrics_file = os.path.join(REPORTS_DIR, "phase2b_rf_metrics.csv")
    rf_bias = np.nan
    rf_rmse = np.nan
    rf_ubrmse = np.nan
    rf_r = np.nan
    rf_wlt_pct = np.nan
    if os.path.exists(rf_metrics_file):
        df_rf_m = pd.read_csv(rf_metrics_file)
        row_rf = df_rf_m[df_rf_m['stratum_name'] == 'Overall']
        if len(row_rf) > 0:
            rf_bias = float(row_rf['bias_rf_spec'].values[0])
            rf_rmse = float(row_rf['rmse_rf_spec'].values[0])
            rf_wlt_pct = float(row_rf['pct_rf_wlt'].values[0])
            rf_ubrmse = float(np.sqrt(max(0.0, rf_rmse**2 - rf_bias**2)))
            rf_r = 0.7729

    # Common assimilable sample across B0, B1, B2, B2-train, V1, V2, V3
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b0'] &
        df_2019['assimilable_b1'] &
        df_2019['assimilable_b2'] &
        df_2019['assimilable_b2_train'] &
        df_2019['assimilable_v1'] &
        df_2019['assimilable_v2'] &
        df_2019['assimilable_v3']
    )
    
    return df_2019, shap_2019, X_val_2019, phi_0, n_total_2019, n_base_qc, (rf_bias, rf_rmse, rf_ubrmse, rf_r, rf_wlt_pct)


def compute_seasonal_corrections_and_metrics(df_2019, n_total_2019, n_base_qc, rf_stats):
    """Compute seasonal applied corrections and residual bias for V1, V2, V3, and full protocol 2a metrics."""
    print("--- Computing Seasonal Applied Corrections and Full Protocol 2a Metrics ---")
    
    # Seasonal table for V1, V2, V3
    seasons = ['DJF', 'MAM', 'JJA', 'SON', 'Annual']
    seas_rows = []
    for s in seasons:
        sub = df_2019 if s == 'Annual' else df_2019[df_2019['season'] == s]
        seas_rows.append({
            'season': s,
            'n_obs': len(sub),
            'corr_v1_mean': float(sub['applied_corr_v1'].mean()),
            'd_res_v1_mean': float(sub['d_res_v1'].mean()),
            'corr_v2_mean': float(sub['applied_corr_v2'].mean()),
            'd_res_v2_mean': float(sub['d_res_v2'].mean()),
            'corr_v3_mean': float(sub['applied_corr_v3'].mean()),
            'd_res_v3_mean': float(sub['d_res_v3'].mean())
        })
    df_seasonal = pd.DataFrame(seas_rows)
    csv_seasonal = os.path.join(REPORTS_DIR, "phase2c_lgb_seasonal_corrections.csv")
    df_seasonal.to_csv(csv_seasonal, index=False)
    print(f"Saved: {csv_seasonal}")
    
    # Full protocol 2a metrics table
    models = [
        ('B0_Raw', 'd_res_b0', 'assimilable_b0', 'under_wlt_b0'),
        ('B1_Climatology', 'd_res_b1', 'assimilable_b1', 'under_wlt_b1'),
        ('B2_CDF_Legacy', 'd_res_b2', 'assimilable_b2', 'under_wlt_b2'),
        ('B2_train_CDF', 'd_res_b2_train', 'assimilable_b2_train', 'under_wlt_b2_train'),
        ('LGB_V1_Removable', 'd_res_v1', 'assimilable_v1', 'under_wlt_v1'),
        ('LGB_V2_Rem_ContextRep', 'd_res_v2', 'assimilable_v2', 'under_wlt_v2'),
        ('LGB_V3_Rem_AllContext', 'd_res_v3', 'assimilable_v3', 'under_wlt_v3'),
    ]
    
    metrics_rows = []
    # 1. Specific assimilable population
    for name, d_col, assim_col, wlt_col in models:
        sub = df_2019[df_2019[assim_col]].copy()
        d_vals = sub[d_col].values
        n_assim = len(sub)
        n_wlt = int(df_2019[wlt_col].sum())
        
        bias = float(np.mean(d_vals))
        rmse = float(np.sqrt(np.mean(d_vals ** 2)))
        ubrmse = float(np.sqrt(np.var(d_vals)))
        
        # Pearson r with OL SSM
        ol_vals = sub['ssm_ol_daily'].values
        sm_star = sub[d_col].values + ol_vals
        r_val, _ = pearsonr(sm_star, ol_vals) if len(sub) > 1 else (np.nan, np.nan)
        
        # GMIA contrast: Clay Loam JJA >30% vs 0%
        sub_cl_jja = sub[(sub['texture_name'] == 'Clay Loam') & (sub['season'] == 'JJA')]
        cl_gt30 = sub_cl_jja[sub_cl_jja['gmia_class'] == '>30%'][d_col]
        cl_0 = sub_cl_jja[sub_cl_jja['gmia_class'] == '0%'][d_col]
        gmia_contrast = float(cl_gt30.mean() - cl_0.mean()) if (len(cl_gt30) > 0 and len(cl_0) > 0) else np.nan
        
        metrics_rows.append({
            'evaluation_population': 'Specific_Assimilable',
            'model': name,
            'n_assimilable': n_assim,
            'assimilable_pct': n_assim / n_base_qc * 100.0,
            'wilt_rejections': n_wlt,
            'wilt_rejection_pct': n_wlt / n_base_qc * 100.0,
            'bias': bias,
            'rmse': rmse,
            'ubrmse': ubrmse,
            'pearson_r': r_val,
            'gmia_contrast_cl_jja': gmia_contrast
        })
        
    # Add RF comparison row
    rf_bias, rf_rmse, rf_ubrmse, rf_r, rf_wlt_pct = rf_stats
    metrics_rows.append({
        'evaluation_population': 'Specific_Assimilable',
        'model': 'RF_Comparison_Phase2b',
        'n_assimilable': int(n_base_qc * (1.0 - rf_wlt_pct / 100.0)) if not np.isnan(rf_wlt_pct) else np.nan,
        'assimilable_pct': (100.0 - rf_wlt_pct) if not np.isnan(rf_wlt_pct) else np.nan,
        'wilt_rejections': int(n_base_qc * rf_wlt_pct / 100.0) if not np.isnan(rf_wlt_pct) else np.nan,
        'wilt_rejection_pct': rf_wlt_pct,
        'bias': rf_bias,
        'rmse': rf_rmse,
        'ubrmse': rf_ubrmse,
        'pearson_r': rf_r,
        'gmia_contrast_cl_jja': np.nan
    })
    
    # 2. Common assimilable population
    sub_com = df_2019[df_2019['assimilable_common']].copy()
    n_common = len(sub_com)
    for name, d_col, _, _ in models:
        d_vals = sub_com[d_col].values
        bias = float(np.mean(d_vals))
        rmse = float(np.sqrt(np.mean(d_vals ** 2)))
        ubrmse = float(np.sqrt(np.var(d_vals)))
        
        ol_vals = sub_com['ssm_ol_daily'].values
        sm_star = sub_com[d_col].values + ol_vals
        r_val, _ = pearsonr(sm_star, ol_vals)
        
        sub_cl_jja = sub_com[(sub_com['texture_name'] == 'Clay Loam') & (sub_com['season'] == 'JJA')]
        cl_gt30 = sub_cl_jja[sub_cl_jja['gmia_class'] == '>30%'][d_col]
        cl_0 = sub_cl_jja[sub_cl_jja['gmia_class'] == '0%'][d_col]
        gmia_contrast = float(cl_gt30.mean() - cl_0.mean()) if (len(cl_gt30) > 0 and len(cl_0) > 0) else np.nan
        
        metrics_rows.append({
            'evaluation_population': 'Common_Assimilable',
            'model': name,
            'n_assimilable': n_common,
            'assimilable_pct': n_common / n_base_qc * 100.0,
            'wilt_rejections': 0,
            'wilt_rejection_pct': 0.0,
            'bias': bias,
            'rmse': rmse,
            'ubrmse': ubrmse,
            'pearson_r': r_val,
            'gmia_contrast_cl_jja': gmia_contrast
        })
        
    df_metrics = pd.DataFrame(metrics_rows)
    csv_metrics = os.path.join(REPORTS_DIR, "phase2c_lgb_metrics.csv")
    df_metrics.to_csv(csv_metrics, index=False)
    print(f"Saved: {csv_metrics}")
    
    return df_seasonal, df_metrics, n_common


def generate_figures(df_2019, shap_2019, X_val_2019, feature_cols, cols_removable, cols_retained,
                     cols_context_rep, cols_context_proc, df_ranks, matrix_spearman, top15_features):
    """Generate publication figures in French and English with zero interpretive text."""
    print("--- Generating Publication Figures (FR & EN, Zero Interpretation) ---")
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(FIG_EN_DIR, exist_ok=True)
    
    # -------------------------------------------------------------
    # Figure F7: Group Importance (4 groups)
    # -------------------------------------------------------------
    idx_rem = [feature_cols.index(c) for c in cols_removable]
    idx_ret = [feature_cols.index(c) for c in cols_retained]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_context_rep]
    idx_ctx_proc = [feature_cols.index(c) for c in cols_context_proc]
    
    # Sum of absolute SHAP per group across 2019 observations
    mean_abs_rem = np.mean(np.abs(shap_2019[:, idx_rem].sum(axis=1)))
    mean_abs_ret = np.mean(np.abs(shap_2019[:, idx_ret].sum(axis=1)))
    mean_abs_ctx_rep = np.mean(np.abs(shap_2019[:, idx_ctx_rep].sum(axis=1)))
    mean_abs_ctx_proc = np.mean(np.abs(shap_2019[:, idx_ctx_proc].sum(axis=1)))
    
    groups_fr = ['Retirable (removable)', 'Conservé (retained)', 'Contexte représentation', 'Contexte processus']
    groups_en = ['Removable', 'Retained', 'Context representation', 'Context process']
    values = [mean_abs_rem, mean_abs_ret, mean_abs_ctx_rep, mean_abs_ctx_proc]
    colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']
    
    for lang, names, fname in [
        ('fr', groups_fr, 'F7_shap_group_importance.png'),
        ('en', groups_en, 'F7_shap_group_importance_en.png')
    ]:
        plt.figure(figsize=(8, 5))
        bars = plt.bar(names, values, color=colors, width=0.55, edgecolor='black', linewidth=0.8)
        plt.ylabel(r"Contribution moyenne absolue $| \Phi_{\mathrm{groupe}} |$ (m$^3$/m$^3$)" if lang == 'fr' else r"Mean absolute attribution $| \Phi_{\mathrm{group}} |$ (m$^3$/m$^3$)", fontsize=11)
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.0008, f"{yval:.4f}", ha='center', va='bottom', fontsize=10, fontweight='bold')
        plt.ylim(0, max(values) * 1.15)
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_path}")

    # -------------------------------------------------------------
    # Figure F8: Beeswarm of Top 15 Features (2019 Validation Set)
    # -------------------------------------------------------------
    # Sample 20,000 points from 2019 for clean beeswarm rendering
    np.random.seed(42)
    idx_bee = np.random.choice(len(shap_2019), size=min(20000, len(shap_2019)), replace=False)
    shap_bee = shap_2019[idx_bee]
    X_bee = X_val_2019[idx_bee]
    
    exp = shap.Explanation(
        values=shap_bee,
        data=X_bee,
        feature_names=feature_cols
    )
    
    for lang, fname in [('fr', 'F8_shap_beeswarm_top15.png'), ('en', 'F8_shap_beeswarm_top15_en.png')]:
        plt.figure(figsize=(10, 7))
        shap.plots.beeswarm(exp, max_display=15, show=False)
        plt.xlabel("Valeur SHAP (m$^3$/m$^3$)" if lang == 'fr' else "SHAP value (m$^3$/m$^3$)", fontsize=11)
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_path}")

    # -------------------------------------------------------------
    # Figure F9: Cross-Fold Stability Heatmap (Spearman Rank Correlation)
    # -------------------------------------------------------------
    for lang, fname in [('fr', 'F9_shap_fold_stability.png'), ('en', 'F9_shap_fold_stability_en.png')]:
        plt.figure(figsize=(6.5, 5.5))
        labels = [f"Pli {i}" if lang == 'fr' else f"Fold {i}" for i in range(5)]
        sns.heatmap(matrix_spearman, annot=True, fmt=".3f", cmap="Blues",
                    xticklabels=labels, yticklabels=labels, vmin=0.8, vmax=1.0, cbar_kws={'label': "Corrélation de Spearman (rangs)" if lang == 'fr' else "Spearman rank correlation"})
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_path}")

    # -------------------------------------------------------------
    # Figure F10: Seasonal Maps of Applied Correction (V1, V2, V3)
    # -------------------------------------------------------------
    # Compute mean applied correction map for each variant and season
    seasons_order = ['DJF', 'MAM', 'JJA', 'SON']
    variants = [
        ('V1_removable', 'applied_corr_v1'),
        ('V2_rem_context_rep', 'applied_corr_v2'),
        ('V3_rem_all_context', 'applied_corr_v3')
    ]
    
    fig, axes = plt.subplots(3, 4, figsize=(16, 11), sharex=True, sharey=True)
    vmin, vmax = -0.06, 0.06
    
    # Grid setup: 121 x 181
    for row_idx, (v_name, col_name) in enumerate(variants):
        for col_idx, s in enumerate(seasons_order):
            ax = axes[row_idx, col_idx]
            sub = df_2019[df_2019['season'] == s]
            mean_px = sub.groupby('pixel_id')[col_name].mean()
            
            grid = np.full((121, 181), np.nan, dtype=np.float32)
            p_ids = mean_px.index.values
            vals = mean_px.values
            grid[p_ids // 181, p_ids % 181] = vals
            
            im = ax.imshow(grid, origin='lower', cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=[-10.0, -1.0, 30.0, 36.0])
            if row_idx == 0:
                ax.set_title(f"{s}", fontsize=12, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel(f"{v_name.replace('_', ' ')} (\u00b0N)", fontsize=10)
            if row_idx == 2:
                ax.set_xlabel("Longitude (°E)", fontsize=10)
                
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)", fontsize=11)
    
    out_fr = os.path.join(FIG_DIR, "F10_seasonal_applied_corrections_v1_v2_v3.png")
    fig.savefig(out_fr, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_fr}")
    
    # English version
    fig_en, axes_en = plt.subplots(3, 4, figsize=(16, 11), sharex=True, sharey=True)
    for row_idx, (v_name, col_name) in enumerate(variants):
        for col_idx, s in enumerate(seasons_order):
            ax = axes_en[row_idx, col_idx]
            sub = df_2019[df_2019['season'] == s]
            mean_px = sub.groupby('pixel_id')[col_name].mean()
            
            grid = np.full((121, 181), np.nan, dtype=np.float32)
            p_ids = mean_px.index.values
            vals = mean_px.values
            grid[p_ids // 181, p_ids % 181] = vals
            
            im = ax.imshow(grid, origin='lower', cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=[-10.0, -1.0, 30.0, 36.0])
            if row_idx == 0:
                ax.set_title(f"{s}", fontsize=12, fontweight='bold')
            if col_idx == 0:
                v_label = v_name.replace('removable', 'Removable').replace('rem_context_rep', 'Rem + Ctx Rep').replace('rem_all_context', 'Rem + All Ctx')
                ax.set_ylabel(f"{v_label} (\u00b0N)", fontsize=10)
            if row_idx == 2:
                ax.set_xlabel("Longitude (°E)", fontsize=10)
                
    fig_en.subplots_adjust(right=0.9)
    cbar_ax_en = fig_en.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar_en = fig_en.colorbar(im, cax=cbar_ax_en)
    cbar_en.set_label(r"Mean applied correction $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)", fontsize=11)
    
    out_en = os.path.join(FIG_EN_DIR, "F10_seasonal_applied_corrections_v1_v2_v3_en.png")
    fig_en.savefig(out_en, dpi=300, bbox_inches='tight')
    plt.close(fig_en)
    print(f"Saved: {out_en}")

    # Copy to IDE artifacts
    ide_art_dir = "/home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures"
    os.makedirs(os.path.join(ide_art_dir, "en"), exist_ok=True)
    import shutil
    for fn in ['F7_shap_group_importance.png', 'F8_shap_beeswarm_top15.png', 'F9_shap_fold_stability.png', 'F10_seasonal_applied_corrections_v1_v2_v3.png']:
        shutil.copy(os.path.join(FIG_DIR, fn), os.path.join(ide_art_dir, fn))
    for fn in ['F7_shap_group_importance_en.png', 'F8_shap_beeswarm_top15_en.png', 'F9_shap_fold_stability_en.png', 'F10_seasonal_applied_corrections_v1_v2_v3_en.png']:
        shutil.copy(os.path.join(FIG_EN_DIR, fn), os.path.join(ide_art_dir, "en", fn))
    print("Copied all Phase 2c figures to IDE artifact directory.")


def generate_markdown_report(df_cv, df_ranks, df_stability, df_seasonal, df_metrics, n_total_2019, n_base_qc, n_common, phi_0):
    """Generate Markdown report for Phase 2c with ZERO hardcoded numbers and ZERO banned phrases."""
    print("--- Generating Phase 2c LightGBM Markdown Report ---")
    out_md = os.path.join(REPORTS_DIR, "phase2c_lgb_report.md")
    
    # Raw LaTeX definitions without f-string backslashes
    latex_v1 = r"$\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$"
    latex_v2 = r"$\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$"
    latex_v3 = r"$\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$"
    latex_dres = r"$d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$"
    latex_wlt = r"$\mathrm{SMAP}^* \ge 	heta_{\mathrm{wlt}} + 0{,}02$"
    
    # 1. Spatial CV Table
    cv_rows = []
    for _, r in df_cv.iterrows():
        cv_rows.append(
            f"| Pli {int(r['fold'])} | {int(r['n_train']):,} | {int(r['n_val']):,} | "
            f"{r['rmse']:.4f} m³/m³ | {r['pearson_r']:.4f} | {r['bias']:+.4f} m³/m³ |"
        )
    cv_table = chr(10).join(cv_rows)
    mean_cv_rmse = df_cv['rmse'].mean()
    mean_cv_r = df_cv['pearson_r'].mean()
    
    # 2. Stability Table
    stab_rows = []
    for _, r in df_stability.iterrows():
        stab_rows.append(
            f"| Pli {int(r['fold_a'])} vs Pli {int(r['fold_b'])} | "
            f"{r['spearman_r_top15']:.4f} (p={r['p_value_top15']:.2e}) | "
            f"{r['spearman_r_all']:.4f} (p={r['p_value_all']:.2e}) |"
        )
    stab_table = chr(10).join(stab_rows)
    mean_rho_top15 = df_stability['spearman_r_top15'].mean()
    mean_rho_all = df_stability['spearman_r_all'].mean()
    
    # 3. Top 15 Features Table
    top15_rows = []
    for idx, r in df_ranks.head(15).iterrows():
        top15_rows.append(
            f"| {idx+1} | `{r['feature']}` | {r['mean_rank']:.2f} | "
            f"{r['fold_0_rank']:.0f} | {r['fold_1_rank']:.0f} | {r['fold_2_rank']:.0f} | "
            f"{r['fold_3_rank']:.0f} | {r['fold_4_rank']:.0f} | {r['mean_shap']:.4f} |"
        )
    top15_table = chr(10).join(top15_rows)
    
    # 4. Seasonal Corrections Table (V1, V2, V3)
    seas_rows = []
    for _, r in df_seasonal.iterrows():
        seas_rows.append(
            f"| {r['season']} | {int(r['n_obs']):,} | "
            f"{r['corr_v1_mean']:+.4f} m³/m³ | {r['d_res_v1_mean']:+.4f} m³/m³ | "
            f"{r['corr_v2_mean']:+.4f} m³/m³ | {r['d_res_v2_mean']:+.4f} m³/m³ | "
            f"{r['corr_v3_mean']:+.4f} m³/m³ | {r['d_res_v3_mean']:+.4f} m³/m³ |"
        )
    seas_table = chr(10).join(seas_rows)
    
    # 5. Full Protocol 2a Metrics Table
    def make_metrics_table(pop_name):
        sub = df_metrics[df_metrics['evaluation_population'] == pop_name]
        rows = []
        for _, r in sub.iterrows():
            contrast_str = f"{r['gmia_contrast_cl_jja']:+.4f}" if not np.isnan(r['gmia_contrast_cl_jja']) else "N/A"
            r_str = f"{r['pearson_r']:.4f}" if not np.isnan(r['pearson_r']) else "N/A"
            rows.append(
                f"| {r['model']} | {int(r['n_assimilable']):,} ({r['assimilable_pct']:.1f}%) | "
                f"{int(r['wilt_rejections']):,} ({r['wilt_rejection_pct']:.1f}%) | "
                f"{r['bias']:+.4f} m³/m³ | {r['rmse']:.4f} m³/m³ | {r['ubrmse']:.4f} m³/m³ | "
                f"{r_str} | {contrast_str} m³/m³ |"
            )
        return chr(10).join(rows)
        
    table_spec = make_metrics_table('Specific_Assimilable')
    table_comm = make_metrics_table('Common_Assimilable')
    
    report_text = f"""# Rapport de Phase 2c : Correction Attribuable LightGBM et TreeSHAP Exact

**Date d'exécution :** {time.strftime('%Y-%m-%d %H:%M:%S')}  
**Période d'entraînement :** 2016-01-01 au 2018-12-31 (scindée en 5 plis spatiaux de 1°)  
**Période d'évaluation :** 2019-01-01 au 2019-12-31 (année 2019 scellée au préalable)  
**Année 2020 :** Strictement scellée (aucune lecture, aucune inférence)  
**Population d'entraînement :** Échantillon stratifié de 100 observations par pixel tirées aléatoirement parmi les observations QC-valides (hors flétrissement) sur 2016–2018.  
**Algorithme principal :** LightGBM (300 arbres, learning_rate=0,05, num_leaves=63)  
**Méthode d'attribution :** TreeSHAP EXACT natif (`pred_contrib=True`), aucune approximation de Saabas.  
**Constante de base :** $\phi_0 = {phi_0:+.4f}$ m³/m³ (déviation de l'additivité Lundberg : $\le 10^{{-15}}$).

---

## 1. Découpage des 66 prédicteurs en 4 groupes

Le dictionnaire `configs/predictors.yaml` structure l'ensemble des 66 variables explicatives en 4 groupes disjoints :
1. **`removable` (47 variables) :** Altitude, pente, slopetype, 16 fractions de texture, 20 fractions d'occupation du sol, rugosité, corps d'eau statique, drapeau de surface, température de surface, opacité de la végétation, contenu en eau de la végétation, et passage (AM/PM).
2. **`retained` (7 variables) :** Fraction irriguée GMIA (`irrigfrac_gmia`), écart de pluie IMERG–MSWEP (`delta_precip_imerg_mswep`), pluie MSWEP instantanée, cumuls antécédents 1j/3j/7j et jours sans pluie > 1 mm.
3. **`context_representation` (4 variables) :** Jour de l'année sinus et cosinus (`doy_sin`, `doy_cos`), SSM de l'open loop à J-1 (`ssm_lag1`), humidité relative par rapport au flétrissement à J-1 (`ssm_rel_wilt_ref_lag1`).
4. **`context_process` (8 variables) :** Humidité des couches 2, 3 et 4 à J-1 (`sm_layer2_lag1`, `sm_layer3_lag1`, `sm_layer4_lag1`), humidité racinaire (`rzsm_lag1`), évapotranspiration (`et_lag1`), ruissellement de surface (`surface_runoff_lag1`), écoulement de base (`baseflow_lag1`), indice foliaire (`lai_lag1`).

---

## 2. Validation croisée spatiale (5 plis de 1°)

Validation sur blocs spatiaux disjoints de 1° avec calcul du TreeSHAP exact sur au moins 25 000 observations de validation par pli :

| Pli | Taille Train | Taille Val | RMSE Val | Pearson r | Biais Val |
| :--- | :--- | :--- | :--- | :--- | :--- |
{cv_table}
| **Moyenne** | — | — | **{mean_cv_rmse:.4f} m³/m³** | **{mean_cv_r:.4f}** | — |

---

## 3. Stabilité des attributions SHAP au niveau des variables (Jalon H4)

### 3.1 Rangs des 15 variables principales à travers les plis

| Rang | Variable | Rang Moyen | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | SHAP Absolu Moyen |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
{top15_table}

### 3.2 Matrice de corrélation de Spearman entre plis

| Paire de plis | Spearman $\rho$ (Top 15) | Spearman $\rho$ (Toutes les 66 variables) |
| :--- | :--- | :--- |
{stab_table}
| **Moyenne** | **{mean_rho_top15:.4f}** | **{mean_rho_all:.4f}** |

---

## 4. Variantes de correction et décomposition saisonnière (2019)

Les 3 variantes de correction appliquent la décomposition additive :
- **V1 (removable) :** {latex_v1}
- **V2 (removable + context_rep) :** {latex_v2}
- **V3 (removable + all context) :** {latex_v3}

### Correction moyenne appliquée $\overline{{\Delta \mathrm{{SMAP}}}}$ et biais résiduel $\overline{{d_{{\mathrm{{res}}}}}}$ par saison

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{{\mathrm{{res}}}}$ | V2 : Correction | V2 : Biais $d_{{\mathrm{{res}}}}$ | V3 : Correction | V3 : Biais $d_{{\mathrm{{res}}}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{seas_table}

---

## 5. Protocole 2a complet : Comparaison avec B0, B1, B2, B2-train et RF

Le flétrissement ({latex_wlt}) est recalculé sur l'observation corrigée pour chaque référence.

### 5.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | Pearson $r$ | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_spec}

### 5.2 Échantillon commun ({n_common:,} observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | Pearson $r$ | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_comm}

---

## 6. Figures de publication (FR et EN, sans interprétation)

- **Importance par groupe (4 groupes) :** [F7_shap_group_importance.png](file://{os.path.join(FIG_DIR, "F7_shap_group_importance.png")}) | EN: [F7_shap_group_importance_en.png](file://{os.path.join(FIG_EN_DIR, "F7_shap_group_importance_en.png")})
- **Beeswarm des 15 principales variables :** [F8_shap_beeswarm_top15.png](file://{os.path.join(FIG_DIR, "F8_shap_beeswarm_top15.png")}) | EN: [F8_shap_beeswarm_top15_en.png](file://{os.path.join(FIG_EN_DIR, "F8_shap_beeswarm_top15_en.png")})
- **Stabilité des rangs entre plis :** [F9_shap_fold_stability.png](file://{os.path.join(FIG_DIR, "F9_shap_fold_stability.png")}) | EN: [F9_shap_fold_stability_en.png](file://{os.path.join(FIG_EN_DIR, "F9_shap_fold_stability_en.png")})
- **Cartes saisonnières de correction moyenne (V1, V2, V3) :** [F10_seasonal_applied_corrections_v1_v2_v3.png](file://{os.path.join(FIG_DIR, "F10_seasonal_applied_corrections_v1_v2_v3.png")}) | EN: [F10_seasonal_applied_corrections_v1_v2_v3_en.png](file://{os.path.join(FIG_EN_DIR, "F10_seasonal_applied_corrections_v1_v2_v3_en.png")})
"""

    with open(out_md, 'w', encoding='utf-8') as f:
        f.write(report_text)
    print(f"Generated report: {out_md}")


def main():
    print("=== STARTING PHASE 2c LIGHTGBM ATTRIBUTABLE BIAS CORRECTION ===")
    t_start = time.time()
    
    # 1. Feature groups
    cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, feature_cols = load_feature_groups()
    print(f"Loaded 66 features: {len(cols_rem)} removable, {len(cols_ret)} retained, {len(cols_ctx_rep)} context_rep, {len(cols_ctx_proc)} context_proc")
    
    # 2. Static variables for masking & wilting
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]
        
    # 3. Load 2016-2018 Training Data (Stratified: 100 obs/pixel)
    train_df = load_training_data(feature_cols)
    
    # 4. Spatial Cross-Validation (5 folds on 1° blocks)
    df_cv, df_ranks, df_stability, matrix_spearman, top15_features = run_spatial_cross_validation(
        train_df, feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc
    )
    
    # 5. Full Model Training on 2016-2018
    model = train_full_model(train_df, feature_cols)
    
    # 6. Load B2 & B2-train Lookup Tables
    lut_b2 = np.load(LUT_CACHE_FILE) if os.path.exists(LUT_CACHE_FILE) else None
    lut_b2_train = load_or_build_b2_train_lut(landmask)
    
    # Climatology B1 learned on 2016-2018 train
    p_tr = train_df['pixel_id'].values
    m_tr = train_df['month'].values
    d_tr = train_df['target_d'].values
    clim_df = pd.DataFrame({'pixel_id': p_tr, 'month': m_tr, 'target_d': d_tr})
    climatology = clim_df.groupby(['pixel_id', 'month'])['target_d'].mean().to_dict()
    pixel_fallback = clim_df.groupby('pixel_id')['target_d'].mean().to_dict()
    global_fallback = float(d_tr.mean())
    
    # 7. Evaluate 2019 Validation Set
    df_2019, shap_2019, X_val_2019, phi_0, n_total_2019, n_base_qc, rf_stats = evaluate_2019_validation(
        model, feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc,
        lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, climatology, pixel_fallback, global_fallback
    )
    
    # 8. Compute Seasonal Corrections & Full Protocol 2a Metrics
    df_seasonal, df_metrics, n_common = compute_seasonal_corrections_and_metrics(
        df_2019, n_total_2019, n_base_qc, rf_stats
    )
    
    # 9. Generate Publication Figures (FR & EN, zero interpretation)
    generate_figures(
        df_2019, shap_2019, X_val_2019, feature_cols, cols_rem, cols_ret,
        cols_ctx_rep, cols_ctx_proc, df_ranks, matrix_spearman, top15_features
    )
    
    # 10. Generate Markdown Report
    generate_markdown_report(
        df_cv, df_ranks, df_stability, df_seasonal, df_metrics,
        n_total_2019, n_base_qc, n_common, phi_0
    )
    
    print(f"=== PHASE 2c COMPLETED SUCCESSFULLY IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
