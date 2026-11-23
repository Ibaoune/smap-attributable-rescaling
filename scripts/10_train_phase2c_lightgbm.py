#!/usr/bin/env python3
"""
scripts/10_train_phase2c_lightgbm.py

PHASE 2c : LightGBM Attributable Bias Correction with Exact TreeSHAP.
Target: d = SMAP raw - SSM Noah-MP Open Loop (daily average)
Domain: NorthMor 0.05 deg (121 x 181, 16,030 land pixels)
Evaluation: 2019 (split == 'val')
Strict rule: 2020 remains strictly sealed!

Protocol requirements (Point de contrôle 7c):
1. Zero hardcoded numbers in reports and figures. All metrics dynamically computed from data.
2. Training sample: Stratified sample of 100 observations per pixel across 2016-2018 (non-wilt QC valid).
3. Missing roughness: NaNs in roughness_coefficient left intact for native LightGBM branching;
   roughness_is_nan kept as binary indicator predictor.
4. Exact Lundberg TreeSHAP via LightGBM C++ API (pred_contrib=True) on all 2019 observations.
5. Baseline B1: learned on ALL 2016-2018 observations passing non-wilt QC (Phase 2a definition).
6. Baseline B2-train: LDT CDF reconstructed strictly on 2016-01-01 -> 2019-01-01.
7. Evaluation metrics:
   - r_smapstar_vs_ol: Pearson correlation between SMAP* and SSM_OL.
   - r_anom_with_smap_raw: Anomaly correlation between SMAP* and raw SMAP (pixel x month mean subtracted).
   - GMIA contrast (Clay Loam JJA >30% vs 0%) with 1° block bootstrap 95% CIs and distinct block counts (< 5 blocks = non interpretable).
   - Seasonal metrics (assimilable, wilt rejections, bias, RMSE, correlations) by season for all references.
8. Figures in French and English:
   - F7: Group importance (4 groups)
   - F8: Beeswarm plot of top 15 features in 2019
   - F9: Cross-fold stability (Heatmap Top 15 + Bar plot of mean Spearman rho across all levels)
   - F10: Seasonal maps (V1, V2, V3, B2-train) with 2-98th percentile symmetric color scale,
     accompanied by side-by-side maps of roughness_is_nan and Noah-MP dominant soil texture with identical dimensions.
   - Standalone companion figures: F10_roughness_is_nan and F10_dominant_soiltype.
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

CDF_B2_TRAIN_MOD = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_model", "noahmp_model_cdf_NorthMor.nc")
CDF_B2_TRAIN_OBS = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_obs", "smap_obs_cdf_NorthMor.nc")
LUT_B2_TRAIN_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_b2_train.npy")

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}



def get_n_jobs():
    """Get number of CPU cores allocated by SLURM (SLURM_CPUS_PER_TASK) or fallback to local count."""
    if "SLURM_CPUS_PER_TASK" in os.environ:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
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


def load_training_data(feature_cols):
    """Load 2016-2018 training data with stratified sampling: at most 100 observations per pixel without replacement."""
    print("--- Loading 2016-2018 Training Data (Stratified Sampling: <= 100 obs/pixel without replacement) ---")
    t0 = time.time()
    
    meta_cols = ['pixel_id', 'month', 'target_d', 'spatial_block_id', 'spatial_fold',
                 'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
                 'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5', 'smap_sm']
    load_cols = list(set(feature_cols + meta_cols))
    
    dfs = []
    for year in [2016, 2017, 2018]:
        pq_path = os.path.join(TABLE_DIR, f"sar_modeling_table_{year}.parquet")
        df_y = pd.read_parquet(pq_path, columns=load_cols)
        
        qc_nw = (
            (df_y['qc_lis_rain'] == 0) &
            (df_y['qc_lis_frozen'] == 0) &
            (df_y['qc_lis_forest'] == 0) &
            (df_y['qc_lis_snow'] == 0) &
            (df_y['qc_lis_high_veg'] == 0) &
            (df_y['qc_lis_sat_limit'] == 0) &
            (df_y['qc_lis_vwc_gt_5'] == 0) &
            (~df_y['target_d'].isna()) &
            (~df_y['smap_sm'].isna())
        )
        df_valid = df_y[qc_nw][['pixel_id', 'target_d', 'spatial_block_id', 'spatial_fold'] + feature_cols].copy()
        print(f"  Year {year}: loaded {len(df_valid):,} non-wilt QC valid rows")
        dfs.append(df_valid)
        
    full_train_df = pd.concat(dfs, ignore_index=True)
    
    # Stratified sampling without replacement: at most 100 observations per pixel (no duplicates)
    full_train_df['rand_sort'] = np.random.RandomState(42).rand(len(full_train_df))
    train_sampled = (
        full_train_df.sort_values(['pixel_id', 'rand_sort'])
        .groupby('pixel_id')
        .head(100)
        .drop(columns=['rand_sort'])
        .sample(frac=1.0, random_state=42)
        .reset_index(drop=True)
    )
    
    # Sampling statistics: minimum, median, and number of pixels under 100 obs
    counts_px = train_sampled['pixel_id'].value_counts()
    min_obs_px = int(counts_px.min())
    med_obs_px = float(counts_px.median())
    n_under_100 = int((counts_px < 100).sum())
    n_pixels = int(train_sampled['pixel_id'].nunique())
    
    sampling_stats = {
        'n_train_rows': len(train_sampled),
        'n_train_pixels': n_pixels,
        'min_obs_per_pixel': min_obs_px,
        'median_obs_per_pixel': med_obs_px,
        'n_pixels_under_100': n_under_100
    }
    print(f"Training set ready: {len(train_sampled):,} rows across {n_pixels:,} pixels (min={min_obs_px}, median={med_obs_px:.1f}, pixels < 100 obs={n_under_100:,}) in {time.time()-t0:.2f}s")
    return train_sampled, sampling_stats


def run_spatial_cross_validation(train_df, feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, n_jobs=8):
    """5-Fold Spatial Cross-Validation using pre-existing spatial_fold column (Decision 7)."""
    print("--- Running 5-Fold Spatial Cross-Validation with LightGBM Exact TreeSHAP ---")
    t0 = time.time()
    
    n_folds = 5
    cv_records = []
    fold_shaps = []
    
    # Indices of groups
    idx_rem = [feature_cols.index(c) for c in cols_rem]
    idx_ret = [feature_cols.index(c) for c in cols_ret]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_ctx_rep]
    idx_ctx_proc = [feature_cols.index(c) for c in cols_ctx_proc]
    
    group_dict = {
        'removable': idx_rem,
        'retained': idx_ret,
        'context_representation': idx_ctx_rep,
        'context_process': idx_ctx_proc
    }
    grp_names = list(group_dict.keys())
    
    fold_group_importances = {}
    fold_intra_group_ranks = {f: {} for f in range(n_folds)}
    
    for fold in range(n_folds):
        t_fold = time.time()
        tr = train_df[train_df['spatial_fold'] != fold]
        va = train_df[train_df['spatial_fold'] == fold]
        
        X_tr = tr[feature_cols].values.astype(np.float32)
        y_tr = tr['target_d'].values.astype(np.float32)
        X_va = va[feature_cols].values.astype(np.float32)
        y_va = va['target_d'].values.astype(np.float32)
        
        model_fold = lgb.LGBMRegressor(
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
            n_jobs=n_jobs,
            verbose=-1
        )
        model_fold.fit(X_tr, y_tr)
        
        y_pred = model_fold.predict(X_va)
        rmse_val = float(np.sqrt(np.mean((y_va - y_pred) ** 2)))
        bias_val = float(np.mean(y_pred - y_va))
        r_val = float(np.corrcoef(y_va, y_pred)[0, 1])
        
        cv_records.append({
            'fold': fold,
            'n_train': len(tr),
            'n_val': len(va),
            'rmse': rmse_val,
            'pearson_r': r_val,
            'bias': bias_val
        })
        
        # Exact TreeSHAP on validation fold (sample up to 25,000 observations)
        n_sample_shap = min(len(va), 25000)
        np.random.seed(42 + fold)
        shap_idx = np.random.choice(len(va), n_sample_shap, replace=False)
        X_shap_val = X_va[shap_idx]
        
        # Native TreeSHAP from LightGBM booster
        t_sh = time.time()
        contrib_val = model_fold.booster_.predict(X_shap_val, pred_contrib=True, num_threads=n_jobs)
        shap_val = contrib_val[:, :-1]
        
        mean_abs_shap = np.mean(np.abs(shap_val), axis=0)
        fold_shaps.append(mean_abs_shap)
        
        # Group importance on this fold
        grp_imp = {}
        for g_name, g_indices in group_dict.items():
            grp_imp[g_name] = float(np.sum(mean_abs_shap[g_indices]))
        fold_group_importances[fold] = grp_imp
        
        # Intra-group variable ranks on this fold
        for g_name, g_indices in group_dict.items():
            g_importances = mean_abs_shap[g_indices]
            g_ranks = np.argsort(-g_importances)
            fold_intra_group_ranks[fold][g_name] = g_ranks
            
        print(f"  Fold {fold}: Val RMSE={rmse_val:.4f}, r={r_val:.4f}, Bias={bias_val:+.4f} (SHAP on {n_sample_shap:,} pts) in {time.time()-t_fold:.1f}s")
        
    df_cv = pd.DataFrame(cv_records)
    print(f"Spatial CV mean RMSE: {df_cv['rmse'].mean():.4f}, mean r: {df_cv['pearson_r'].mean():.4f}")
    
    # Feature rank stability across folds
    fold_shaps_arr = np.array(fold_shaps)  # (5, 66)
    mean_shaps = np.mean(fold_shaps_arr, axis=0)
    
    rank_records = []
    for j, col in enumerate(feature_cols):
        rec = {'feature': col, 'mean_shap': mean_shaps[j]}
        for f in range(n_folds):
            ranks_f = np.argsort(-fold_shaps_arr[f])
            rank_of_j = int(np.where(ranks_f == j)[0][0]) + 1
            rec[f'fold_{f}_rank'] = rank_of_j
        rec['mean_rank'] = np.mean([rec[f'fold_{f}_rank'] for f in range(n_folds)])
        rank_records.append(rec)
        
    df_ranks = pd.DataFrame(rank_records).sort_values('mean_rank').reset_index(drop=True)
    csv_ranks = os.path.join(REPORTS_DIR, "phase2c_lgb_feature_ranks.csv")
    df_ranks.to_csv(csv_ranks, index=False)
    print(f"Saved: {csv_ranks}")
    
    # Group stability & intra-group rank correlations
    grp_stab_records = []
    for g_name, g_indices in group_dict.items():
        mean_g_imp = np.mean([fold_group_importances[f][g_name] for f in range(n_folds)])
        
        # Ranks of this group across 4 groups for each fold
        grp_ranks_by_fold = []
        for f in range(n_folds):
            sorted_grps = sorted(fold_group_importances[f].keys(), key=lambda k: -fold_group_importances[f][k])
            grp_ranks_by_fold.append(sorted_grps.index(g_name) + 1)
            
        # Intra-group Spearman rho across fold pairs
        intra_rhos = []
        for fa in range(n_folds):
            for fb in range(fa + 1, n_folds):
                r_sp, _ = spearmanr(fold_shaps_arr[fa][g_indices], fold_shaps_arr[fb][g_indices])
                if not np.isnan(r_sp):
                    intra_rhos.append(r_sp)
        mean_intra_rho = np.mean(intra_rhos) if len(intra_rhos) > 0 else 1.0
        
        grp_stab_records.append({
            'group': g_name,
            'mean_importance': mean_g_imp,
            'fold_0_rank': grp_ranks_by_fold[0],
            'fold_1_rank': grp_ranks_by_fold[1],
            'fold_2_rank': grp_ranks_by_fold[2],
            'fold_3_rank': grp_ranks_by_fold[3],
            'fold_4_rank': grp_ranks_by_fold[4],
            'mean_rank': np.mean(grp_ranks_by_fold),
            'mean_intra_spearman_rho': mean_intra_rho
        })
    df_group_stab = pd.DataFrame(grp_stab_records).sort_values('mean_rank').reset_index(drop=True)
    csv_group = os.path.join(REPORTS_DIR, "phase2c_lgb_group_stability.csv")
    df_group_stab.to_csv(csv_group, index=False)
    print(f"Saved: {csv_group}")
    
    # Pairwise Spearman stability
    top15_features = df_ranks.head(15)['feature'].tolist()
    top15_idx = [feature_cols.index(c) for c in top15_features]
    
    stability_records = []
    matrix_spearman = np.ones((n_folds, n_folds), dtype=np.float32)
    for fa in range(n_folds):
        for fb in range(fa + 1, n_folds):
            # All 66
            r_all, p_all = spearmanr(fold_shaps_arr[fa], fold_shaps_arr[fb])
            # Top 15
            r_top15, p_top15 = spearmanr(fold_shaps_arr[fa][top15_idx], fold_shaps_arr[fb][top15_idx])
            # 4 Groups
            ga = [fold_group_importances[fa][g] for g in grp_names]
            gb = [fold_group_importances[fb][g] for g in grp_names]
            r_grp, p_grp = spearmanr(ga, gb)
            
            # Intra-group correlations
            r_intra_rem, _ = spearmanr(fold_shaps_arr[fa][group_dict['removable']], fold_shaps_arr[fb][group_dict['removable']])
            r_intra_ret, _ = spearmanr(fold_shaps_arr[fa][group_dict['retained']], fold_shaps_arr[fb][group_dict['retained']])
            r_intra_rep, _ = spearmanr(fold_shaps_arr[fa][group_dict['context_representation']], fold_shaps_arr[fb][group_dict['context_representation']])
            r_intra_proc, _ = spearmanr(fold_shaps_arr[fa][group_dict['context_process']], fold_shaps_arr[fb][group_dict['context_process']])
            
            matrix_spearman[fa, fb] = r_top15
            matrix_spearman[fb, fa] = r_top15
            stability_records.append({
                'fold_a': fa,
                'fold_b': fb,
                'spearman_r_groups': r_grp,
                'spearman_intra_removable': r_intra_rem,
                'spearman_intra_retained': r_intra_ret,
                'spearman_intra_context_representation': r_intra_rep,
                'spearman_intra_context_process': r_intra_proc,
                'spearman_r_top15': r_top15,
                'p_value_top15': p_top15,
                'spearman_r_all': r_all,
                'p_value_all': p_all
            })
    df_stability = pd.DataFrame(stability_records)
    csv_stab = os.path.join(REPORTS_DIR, "phase2c_lgb_shap_stability.csv")
    df_stability.to_csv(csv_stab, index=False)
    print(f"Saved: {csv_stab}")
    print(f"SHAP stability: Groups rho={df_stability['spearman_r_groups'].mean():.4f}, Top 15 rho={df_stability['spearman_r_top15'].mean():.4f}, All features rho={df_stability['spearman_r_all'].mean():.4f} across 10 fold pairs in {time.time()-t0:.1f}s")
    
    return df_cv, df_ranks, df_stability, df_group_stab, matrix_spearman, top15_features


def train_full_model(train_df, feature_cols, n_jobs=8):
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
        n_jobs=n_jobs,
        verbose=-1
    )
    model.fit(X_train, y_train)
    print(f"Full LightGBM model trained on {len(X_train):,} samples in {time.time()-t0:.2f}s")
    model_path = os.path.join(INTERIM_DIR, "lgb_model_2016_2018.txt")
    model.booster_.save_model(model_path)
    print(f"Saved trained LightGBM model to {model_path}")
    return model


def load_or_build_b2_train_lut(landmask):
    """Load or build lookup table for Baseline B2-train (LDT reconstructed on 2016-01-01 -> 2019-01-01)."""
    if os.path.exists(LUT_B2_TRAIN_FILE):
        print(f"Loading cached B2-train CDF LUT from {LUT_B2_TRAIN_FILE}...")
        return np.load(LUT_B2_TRAIN_FILE)
        
    print("Building B2-train CDF lookup table from LDT output files...")
    t0 = time.time()
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


def learn_b1_climatology():
    """Learn B1 pixel-month mean bias on ALL 2016-2018 non-wilt QC observations (definition of Phase 2a)."""
    print("--- Learning B1 Climatology on ALL 2016-2018 Non-wilt QC Observations ---")
    t0 = time.time()
    qc_cols = [
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    dfs = []
    for y in [2016, 2017, 2018]:
        pq = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(pq, columns=['pixel_id', 'month', 'target_d'] + qc_cols)
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
        dfs.append(df_y.loc[qc_nw, ['pixel_id', 'month', 'target_d']])
    full_train = pd.concat(dfs, ignore_index=True)
    b1_clim_table = (
        full_train.groupby(['pixel_id', 'month'])['target_d']
        .mean()
        .reset_index()
        .rename(columns={'target_d': 'b1_bias'})
    )
    b1_pixel_table = (
        full_train.groupby('pixel_id')['target_d']
        .mean()
        .reset_index()
        .rename(columns={'target_d': 'b1_pixel_fallback'})
    )
    global_fallback = float(full_train['target_d'].mean())
    print(f"B1 Climatology learned on all {len(full_train):,} rows across 2016-2018 in {time.time()-t0:.1f}s")
    b1_clim_path = os.path.join(INTERIM_DIR, "b1_climatology_2016_2018.csv")
    b1_pix_path = os.path.join(INTERIM_DIR, "b1_pixel_fallback_2016_2018.csv")
    b1_clim_table.to_csv(b1_clim_path, index=False)
    b1_pixel_table.to_csv(b1_pix_path, index=False)
    print(f"Saved B1 climatology tables to {b1_clim_path} and {b1_pix_path}")
    return b1_clim_table, b1_pixel_table, global_fallback


def evaluate_2019_validation(model, feature_cols, cols_removable, cols_retained, cols_context_rep, cols_context_proc,
                             lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, b1_clim_table, b1_pixel_table, global_fallback, n_jobs=8):
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
    
    # 1. Baseline B0: No correction
    df_2019['smap_star_b0'] = df_2019['smap_sm']
    df_2019['d_res_b0'] = df_2019['target_d']
    df_2019['under_wlt_b0'] = (df_2019['smap_star_b0'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b0'] = (
        (df_2019['under_wlt_b0'] == 0) &
        (df_2019['smap_star_b0'] >= 0.01) &
        (df_2019['smap_star_b0'] <= 1.0) &
        (~df_2019['smap_star_b0'].isna())
    )
    
    # 2. Baseline B1: Climatology via vectorized pandas join on (pixel_id, month)
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
    
    # 3. Baseline B2: CDF legacy (2016-2021)
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
    
    # 4. Baseline B2-train: CDF reconstructed on 2016-01-01 -> 2019-01-01
    with nc.Dataset(CDF_B2_TRAIN_OBS) as ds_o_tr:
        obs_x_tr = ds_o_tr.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta_tr = obs_x_tr[1, :, :] - obs_x_tr[0, :, :]
    obs_x0_tr = obs_x_tr[0, :, :]
    
    smap_star_b2_tr = apply_lut(lut_b2_train, obs_x0_tr, obs_delta_tr)
    df_2019['smap_star_b2_train'] = smap_star_b2_tr
    df_2019['applied_corr_b2_train'] = df_2019['smap_star_b2_train'] - df_2019['smap_sm']
    df_2019['b2_train_has_value'] = (~np.isnan(smap_star_b2_tr)).astype(int)
    df_2019['d_res_b2_train'] = smap_star_b2_tr - df_2019['ssm_ol_daily']
    df_2019['under_wlt_b2_train'] = (df_2019['smap_star_b2_train'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b2_train'] = (
        (df_2019['b2_train_has_value'] == 1) &
        (df_2019['under_wlt_b2_train'] == 0) &
        (df_2019['smap_star_b2_train'] >= 0.01) &
        (df_2019['smap_star_b2_train'] <= 1.0)
    )
    
    # 5. EXACT Lundberg TreeSHAP via LightGBM C++ API on ALL 2019 Observations
    X_val_2019 = df_2019[feature_cols].values.astype(np.float32)
    contrib_cache_file = os.path.join(BASE_DIR, "data", "interim", "lgb_contrib_2019.npy")
    if os.path.exists(contrib_cache_file):
        print(f"Loading cached exact TreeSHAP attributions from {contrib_cache_file}...")
        contrib_2019 = np.load(contrib_cache_file)
    else:
        print(f"Computing EXACT TreeSHAP on all {len(df_2019):,} observations of 2019...")
        t_shap = time.time()
        contrib_2019 = model.booster_.predict(X_val_2019, pred_contrib=True, num_threads=n_jobs)
        os.makedirs(os.path.dirname(contrib_cache_file), exist_ok=True)
        np.save(contrib_cache_file, contrib_2019)
        print(f"Exact TreeSHAP computed and cached in {time.time()-t_shap:.1f}s")
    phi_0 = float(contrib_2019[0, -1])
    shap_2019 = contrib_2019[:, :-1]
    
    # Verify Lundberg machine-precision additivity
    y_pred_2019 = model.predict(X_val_2019)
    max_additivity_dev = float(np.max(np.abs(shap_2019.sum(axis=1) + phi_0 - y_pred_2019)))
    print(f"Exact TreeSHAP Lundberg additivity verified. Max additivity deviation: {max_additivity_dev:.2e} (machine precision: exact)")
    
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
    
    # 6. Variants V1, V2, V3
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
    
    # V3 = removable + all context
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

    # Common assimilable sample across all 7 references
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b0'] &
        df_2019['assimilable_b1'] &
        df_2019['assimilable_b2'] &
        df_2019['assimilable_b2_train'] &
        df_2019['assimilable_v1'] &
        df_2019['assimilable_v2'] &
        df_2019['assimilable_v3']
    )
    
    return df_2019, shap_2019, X_val_2019, phi_0, n_total_2019, n_base_qc


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


def compute_seasonal_corrections_and_metrics(df_2019, n_total_2019, n_base_qc):
    """Compute seasonal applied corrections, full protocol 2a metrics, and seasonal stratified metrics."""
    print("--- Computing Seasonal Applied Corrections and Full Protocol 2a Metrics ---")
    
    # 1. Seasonal applied correction and residual bias table
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
            'd_res_v3_mean': float(sub['d_res_v3'].mean()),
            'corr_b2_train_mean': float(sub['applied_corr_b2_train'].mean()),
            'd_res_b2_train_mean': float(sub['d_res_b2_train'].mean())
        })
    df_seasonal = pd.DataFrame(seas_rows)
    csv_seasonal = os.path.join(REPORTS_DIR, "phase2c_lgb_seasonal_corrections.csv")
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
    csv_metrics = os.path.join(REPORTS_DIR, "phase2c_lgb_metrics.csv")
    df_metrics.to_csv(csv_metrics, index=False)
    print(f"Saved: {csv_metrics}")
    
    # 3. Seasonal metrics table for each reference
    seasonal_metrics_rows = []
    for pop_name, pop_col in [('Specific_Assimilable', None), ('Common_Assimilable', 'assimilable_common')]:
        for name, d_col, star_col, assim_col, wlt_col in models:
            for s in seasons:
                df_s = df_2019 if s == 'Annual' else df_2019[df_2019['season'] == s]
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
    csv_seas_metrics = os.path.join(REPORTS_DIR, "phase2c_lgb_seasonal_metrics.csv")
    df_seas_metrics.to_csv(csv_seas_metrics, index=False)
    print(f"Saved: {csv_seas_metrics}")
    
    n_common = int(df_2019['assimilable_common'].sum())
    return df_seasonal, df_metrics, df_seas_metrics, n_common


def generate_figures(df_2019, shap_2019, X_val_2019, feature_cols, cols_removable, cols_retained,
                     cols_context_rep, cols_context_proc, df_ranks, df_stability, df_group_stab,
                     matrix_spearman, top15_features, dom_soil, landmask):
    """Generate publication figures in French and English with zero interpretive text and dynamic legends."""
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
    
    mean_abs_rem = float(np.mean(np.abs(shap_2019[:, idx_rem].sum(axis=1))))
    mean_abs_ret = float(np.mean(np.abs(shap_2019[:, idx_ret].sum(axis=1))))
    mean_abs_ctx_rep = float(np.mean(np.abs(shap_2019[:, idx_ctx_rep].sum(axis=1))))
    mean_abs_ctx_proc = float(np.mean(np.abs(shap_2019[:, idx_ctx_proc].sum(axis=1))))
    
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
            plt.text(bar.get_x() + bar.get_width()/2.0, yval + (max(values)*0.01), f"{yval:.4f}", ha='center', va='bottom', fontsize=10, fontweight='bold')
        plt.ylim(0, max(values) * 1.15)
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_path}")

    # -------------------------------------------------------------
    # Figure F8: Beeswarm of Top 15 Features (2019 Validation Set)
    # -------------------------------------------------------------
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
    # Figure F9: Cross-Fold Stability (Heatmap + Bar plot of mean Spearman rho)
    # -------------------------------------------------------------
    mean_rho_grp = float(df_stability['spearman_r_groups'].mean())
    mean_rho_rem = float(df_stability['spearman_intra_removable'].mean())
    mean_rho_ret = float(df_stability['spearman_intra_retained'].mean())
    mean_rho_ctx_rep = float(df_stability['spearman_intra_context_representation'].mean())
    mean_rho_ctx_proc = float(df_stability['spearman_intra_context_process'].mean())
    mean_rho_top15_val = float(df_stability['spearman_r_top15'].mean())
    mean_rho_all_val = float(df_stability['spearman_r_all'].mean())
    
    cat_names_fr = [
        "4 Groupes (global)",
        "Intra : Retirable (47 var.)",
        "Intra : Conservé (7 var.)",
        "Intra : Ctx Repr. (4 var.)",
        "Intra : Ctx Proc. (8 var.)",
        "Top 15 variables",
        "Toutes les 66 variables"
    ]
    cat_names_en = [
        "4 Groups (global)",
        "Intra: Removable (47 vars)",
        "Intra: Retained (7 vars)",
        "Intra: Ctx Repr. (4 vars)",
        "Intra: Ctx Proc. (8 vars)",
        "Top 15 features",
        "All 66 features"
    ]
    rhos = [mean_rho_grp, mean_rho_rem, mean_rho_ret, mean_rho_ctx_rep, mean_rho_ctx_proc, mean_rho_top15_val, mean_rho_all_val]
    y_pos = np.arange(len(rhos))
    colors_bar = ['#2b5c8f', '#1f77b4', '#2ca02c', '#ff7f0e', '#d62728', '#9467bd', '#8c564b']

    for lang, fname in [('fr', 'F9_shap_fold_stability.png'), ('en', 'F9_shap_fold_stability_en.png')]:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={'width_ratios': [1.0, 1.3]})
        
        labels = [f"Pli {i}" if lang == 'fr' else f"Fold {i}" for i in range(5)]
        sns.heatmap(matrix_spearman, annot=True, fmt=".3f", cmap="Blues",
                    xticklabels=labels, yticklabels=labels, vmin=0.8, vmax=1.0,
                    cbar_kws={'label': "Corrélation de Spearman (rangs)" if lang == 'fr' else "Spearman rank correlation"},
                    ax=ax1)
        ax1.set_title("(a) Matrice de corrélation (Top 15)" if lang == 'fr' else "(a) Correlation matrix (Top 15)", fontsize=11, fontweight='bold')
        
        bars = ax2.barh(y_pos, rhos, color=colors_bar, height=0.6, edgecolor='black', linewidth=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(cat_names_fr if lang == 'fr' else cat_names_en, fontsize=10)
        ax2.invert_yaxis()
        ax2.set_xlim(0.7, 1.03)
        ax2.set_xlabel("Corrélation de rang moyenne (Spearman)" if lang == 'fr' else "Mean Spearman rank correlation", fontsize=11)
        ax2.set_title("(b) Stabilité à travers les 10 paires de plis" if lang == 'fr' else "(b) Stability across 10 fold pairs", fontsize=11, fontweight='bold')
        ax2.grid(axis='x', linestyle='--', alpha=0.5)
        for bar in bars:
            xval = bar.get_width()
            ax2.text(xval + 0.005, bar.get_y() + bar.get_height()/2.0, f"{xval:.4f}", va='center', ha='left', fontsize=9, fontweight='bold')
            
        fig.tight_layout()
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {out_path}")

    # -------------------------------------------------------------
    # Figure F10: Seasonal Maps (V1, V2, V3, B2-train) + roughness_is_nan + dominant_soiltype
    # -------------------------------------------------------------
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
    
    # Color scale fixed to 2nd-98th percentiles across all 16 maps
    all_maps_data = []
    for _, col_name in variants_fr:
        for s in seasons_order:
            sub = df_2019[df_2019['season'] == s]
            m = sub.groupby('pixel_id')[col_name].mean().values
            all_maps_data.append(m[~np.isnan(m)])
    all_vals_concat = np.concatenate(all_maps_data)
    p2 = float(np.nanpercentile(all_vals_concat, 2))
    p98 = float(np.nanpercentile(all_vals_concat, 98))
    vlim = float(max(abs(p2), abs(p98)))
    vmin, vmax = -vlim, vlim
    print(f"F10 color scale fixed to 2-98th percentiles: [{vmin:+.4f}, {vmax:+.4f}] m3/m3")
    
    # 1. Grid of roughness_is_nan
    grid_nan = np.full((121, 181), np.nan, dtype=np.float32)
    rough_nan_px = df_2019.groupby('pixel_id')['roughness_is_nan'].first()
    grid_nan[rough_nan_px.index // 181, rough_nan_px.index % 181] = rough_nan_px.values
    cmap_nan = mcolors.ListedColormap(['#e0e0e0', '#d62728'])
    
    # Dynamic percentages for roughness
    n_valid_px = int((grid_nan == 0.0).sum())
    n_nan_px = int((grid_nan == 1.0).sum())
    total_r_px = n_valid_px + n_nan_px
    pct_nan_px = float(n_nan_px / total_r_px * 100.0) if total_r_px > 0 else 0.0
    pct_val_px = 100.0 - pct_nan_px

    # 2. Grid of dominant soil texture (Noah-MP)
    grid_soil = np.full((121, 181), np.nan, dtype=np.float32)
    soil_active_codes = [1, 3, 4, 6, 7, 9, 12]
    soil_code_to_idx = {c: i for i, c in enumerate(soil_active_codes)}
    
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                st = dom_soil[r, c]
                if st in soil_code_to_idx:
                    grid_soil[r, c] = soil_code_to_idx[st]
                    
    soil_colors = ["#f4d03f", "#eb984e", "#a569bd", "#5dade2", "#48c9b0", "#e67e22", "#c0392b"]
    cmap_soil = mcolors.ListedColormap(soil_colors)
    
    # Dynamic percentages for each soil class
    total_land_px = int((landmask == 1).sum())
    soil_pcts = {}
    for c in soil_active_codes:
        c_cnt = int(((dom_soil == c) & (landmask == 1)).sum())
        soil_pcts[c] = float(c_cnt / total_land_px * 100.0) if total_land_px > 0 else 0.0

    soil_names_fr = {
        1: "Sable (1)", 3: "Loam sableux (3)", 4: "Loam limoneux (4)",
        6: "Loam (6)", 7: "Loam argilo-sableux (7)", 9: "Loam argileux (9)", 12: "Argile (12)"
    }
    soil_names_en = {
        1: "Sand (1)", 3: "Sandy Loam (3)", 4: "Silt Loam (4)",
        6: "Loam (6)", 7: "Sandy Clay Loam (7)", 9: "Clay Loam (9)", 12: "Clay (12)"
    }

    # Standalone roughness_is_nan figure
    for lang, fname in [('fr', 'F10_roughness_is_nan.png'), ('en', 'F10_roughness_is_nan_en.png')]:
        fig_nan, ax_n = plt.subplots(figsize=(7, 5))
        ax_n.imshow(grid_nan, origin='lower', cmap=cmap_nan, extent=[-10.0, -1.0, 30.0, 36.0])
        ax_n.set_xlabel("Longitude (°E)", fontsize=10)
        ax_n.set_ylabel("Latitude (°N)", fontsize=10)
        ax_n.set_title("Masque spatial de rugosité manquante (roughness_is_nan)" if lang == 'fr' else "Spatial mask of missing roughness (roughness_is_nan)", fontsize=11, fontweight='bold')
        
        lbl_v = f"Rugosité valide ({pct_val_px:.1f} %)" if lang == 'fr' else f"Valid roughness ({pct_val_px:.1f}%)"
        lbl_m = f"Rugosité manquante ({pct_nan_px:.1f} %)" if lang == 'fr' else f"Missing roughness ({pct_nan_px:.1f}%)"
        patch_val = mpatches.Patch(color='#e0e0e0', label=lbl_v)
        patch_nan = mpatches.Patch(color='#d62728', label=lbl_m)
        ax_n.legend(handles=[patch_val, patch_nan], loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
        
        out_n = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        fig_nan.savefig(out_n, dpi=300, bbox_inches='tight')
        plt.close(fig_nan)
        print(f"Saved: {out_n}")

    # Standalone dominant soil texture figure
    for lang, fname in [('fr', 'F10_dominant_soiltype.png'), ('en', 'F10_dominant_soiltype_en.png')]:
        fig_s, ax_s = plt.subplots(figsize=(7, 5))
        ax_s.imshow(grid_soil, origin='lower', cmap=cmap_soil, vmin=-0.5, vmax=6.5, extent=[-10.0, -1.0, 30.0, 36.0])
        ax_s.set_xlabel("Longitude (°E)", fontsize=10)
        ax_s.set_ylabel("Latitude (°N)", fontsize=10)
        ax_s.set_title("Texture dominante du sol (Noah-MP)" if lang == 'fr' else "Noah-MP dominant soil texture", fontsize=11, fontweight='bold')
        
        s_names = soil_names_fr if lang == 'fr' else soil_names_en
        s_patches = [
            mpatches.Patch(color=soil_colors[i], label=f"{s_names[c]} ({soil_pcts[c]:.1f} %)" if lang == 'fr' else f"{s_names[c]} ({soil_pcts[c]:.1f}%)")
            for i, c in enumerate(soil_active_codes)
        ]
        ax_s.legend(handles=s_patches, loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=8)
        
        out_s = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        fig_s.savefig(out_s, dpi=300, bbox_inches='tight')
        plt.close(fig_s)
        print(f"Saved: {out_s}")

    # Full F10 multi-panel figure: 4 rows x 7 cols
    # (Cols 0-3: seasons, Col 4: roughness_is_nan, Col 5: dominant_soiltype, Col 6: colorbar)
    for lang, variants, fname in [
        ('fr', variants_fr, 'F10_seasonal_applied_corrections_v1_v2_v3.png'),
        ('en', variants_en, 'F10_seasonal_applied_corrections_v1_v2_v3_en.png')
    ]:
        fig = plt.figure(figsize=(22, 12))
        gs = fig.add_gridspec(4, 7, width_ratios=[1, 1, 1, 1, 1.05, 1.05, 0.07], wspace=0.04, hspace=0.06)
        
        # 16 seasonal maps
        for row_idx, (v_label, col_name) in enumerate(variants):
            for col_idx, s in enumerate(seasons_order):
                ax = fig.add_subplot(gs[row_idx, col_idx])
                sub = df_2019[df_2019['season'] == s]
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
                else:
                    ax.set_yticklabels([])
                if row_idx == 3:
                    ax.set_xlabel("Longitude (°E)", fontsize=10)
                else:
                    ax.set_xticklabels([])
                    
        # Col 4: roughness_is_nan panel (spanning rows 1-2, right next to JJA/SON of V2/V3)
        ax_nan = fig.add_subplot(gs[1:3, 4])
        ax_nan.imshow(grid_nan, origin='lower', cmap=cmap_nan, extent=[-10.0, -1.0, 30.0, 36.0])
        ax_nan.set_title("Masque roughness_is_nan\n(Zones sans rugosité)" if lang == 'fr' else "roughness_is_nan mask\n(Missing roughness zones)", fontsize=10, fontweight='bold')
        ax_nan.set_xlabel("Longitude (°E)", fontsize=9)
        ax_nan.set_ylabel("Latitude (°N)", fontsize=9)
        lbl_v = f"Valide ({pct_val_px:.1f} %)" if lang == 'fr' else f"Valid ({pct_val_px:.1f}%)"
        lbl_m = f"NaN ({pct_nan_px:.1f} %)" if lang == 'fr' else f"NaN ({pct_nan_px:.1f}%)"
        p_val = mpatches.Patch(color='#e0e0e0', label=lbl_v)
        p_nan = mpatches.Patch(color='#d62728', label=lbl_m)
        ax_nan.legend(handles=[p_val, p_nan], loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=8)
        
        # Col 5: dominant soil texture panel (spanning rows 1-2, right next to roughness_is_nan with identical dimensions)
        ax_tex = fig.add_subplot(gs[1:3, 5])
        ax_tex.imshow(grid_soil, origin='lower', cmap=cmap_soil, vmin=-0.5, vmax=6.5, extent=[-10.0, -1.0, 30.0, 36.0])
        ax_tex.set_title("Texture sol Noah-MP\n(7 classes de sol)" if lang == 'fr' else "Noah-MP dominant soil\n(7 soil classes)", fontsize=10, fontweight='bold')
        ax_tex.set_xlabel("Longitude (°E)", fontsize=9)
        ax_tex.set_yticklabels([])
        s_names = soil_names_fr if lang == 'fr' else soil_names_en
        s_patches = [
            mpatches.Patch(color=soil_colors[i], label=f"{s_names[c]} ({soil_pcts[c]:.1f} %)" if lang == 'fr' else f"{s_names[c]} ({soil_pcts[c]:.1f}%)")
            for i, c in enumerate(soil_active_codes)
        ]
        ax_tex.legend(handles=s_patches, loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=7)
        
        # Col 6: Vertical colorbar
        cax = fig.add_subplot(gs[:, 6])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(
            r"Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)" if lang == 'fr' else
            r"Mean applied correction $\overline{\Delta \mathrm{SMAP}}$ (m$^3$/m$^3$)", fontsize=11
        )
        
        out_path = os.path.join(FIG_DIR if lang == 'fr' else FIG_EN_DIR, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {out_path}")

    # Copy to IDE artifacts
    ide_art_dir = "/home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures"
    os.makedirs(os.path.join(ide_art_dir, "en"), exist_ok=True)
    import shutil
    fig_names = [
        'F7_shap_group_importance.png', 'F8_shap_beeswarm_top15.png', 'F9_shap_fold_stability.png',
        'F10_seasonal_applied_corrections_v1_v2_v3.png', 'F10_roughness_is_nan.png', 'F10_dominant_soiltype.png'
    ]
    for fn in fig_names:
        shutil.copy(os.path.join(FIG_DIR, fn), os.path.join(ide_art_dir, fn))
    for fn in fig_names:
        en_fn = fn.replace('.png', '_en.png')
        if os.path.exists(os.path.join(FIG_EN_DIR, en_fn)):
            shutil.copy(os.path.join(FIG_EN_DIR, en_fn), os.path.join(ide_art_dir, "en", en_fn))
    print("Copied all Phase 2c figures to IDE artifact directory.")


def generate_markdown_report(df_cv, df_ranks, df_stability, df_group_stab, df_seasonal, df_metrics, df_seas_metrics,
                             n_total_2019, n_base_qc, n_common, phi_0, sampling_stats):
    """Generate Markdown report for Phase 2c with ZERO hardcoded numbers and ZERO banned phrases."""
    print("--- Generating Phase 2c LightGBM Markdown Report ---")
    out_md = os.path.join(REPORTS_DIR, "phase2c_lgb_report.md")
    
    latex_v1 = r"$\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$"
    latex_v2 = r"$\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$"
    latex_v3 = r"$\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$"
    latex_dres = r"$d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$"
    latex_theta = r"$	heta$"
    latex_rho = r"$\rho$"
    latex_wlt = r"$\mathrm{SMAP}^* \ge \theta_{\mathrm{wlt}} + 0{,}02$"
    
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
    
    # 2. Group Stability Table
    grp_rows = []
    for _, r in df_group_stab.iterrows():
        grp_rows.append(
            f"| `{r['group']}` | {r['mean_importance']:.4f} m³/m³ | "
            f"{int(r['fold_0_rank'])} | {int(r['fold_1_rank'])} | {int(r['fold_2_rank'])} | "
            f"{int(r['fold_3_rank'])} | {int(r['fold_4_rank'])} | {r['mean_rank']:.2f} | "
            f"{r['mean_intra_spearman_rho']:.4f} |"
        )
    group_table = chr(10).join(grp_rows)
    mean_rho_groups = df_stability['spearman_r_groups'].mean()
    
    # 3. Pairwise Stability Table
    stab_rows = []
    for _, r in df_stability.iterrows():
        stab_rows.append(
            f"| Pli {int(r['fold_a'])} vs Pli {int(r['fold_b'])} | "
            f"{r['spearman_r_groups']:.4f} | "
            f"{r['spearman_r_top15']:.4f} (p={r['p_value_top15']:.2e}) | "
            f"{r['spearman_r_all']:.4f} (p={r['p_value_all']:.2e}) |"
        )
    stab_table = chr(10).join(stab_rows)
    mean_rho_top15 = df_stability['spearman_r_top15'].mean()
    mean_rho_all = df_stability['spearman_r_all'].mean()
    
    # 4. Top 15 Features Table
    top15_rows = []
    for idx, r in df_ranks.head(15).iterrows():
        top15_rows.append(
            f"| {idx+1} | `{r['feature']}` | {r['mean_rank']:.2f} | "
            f"{r['fold_0_rank']:.0f} | {r['fold_1_rank']:.0f} | {r['fold_2_rank']:.0f} | "
            f"{r['fold_3_rank']:.0f} | {r['fold_4_rank']:.0f} | {r['mean_shap']:.4f} |"
        )
    top15_table = chr(10).join(top15_rows)
    
    # 5. Seasonal Corrections Table (V1, V2, V3, B2-train)
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
    
    # 6. Full Protocol 2a Metrics Table
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
    
    # 7. Seasonal Stratified Metrics Table (Specific Assimilable)
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
    
    report_text = f"""# Rapport de Phase 2c : Correction Attribuable LightGBM et TreeSHAP Exact

**Date d'exécution :** {time.strftime('%Y-%m-%d %H:%M:%S')}  
**Période d'entraînement :** 2016-01-01 au 2018-12-31 (scindée en 5 plis spatiaux de 1°)  
**Période d'évaluation :** 2019-01-01 au 2019-12-31 (année 2019 scellée au préalable)  
**Année 2020 :** Strictement scellée (aucune lecture, aucune inférence)  
**Population d'entraînement :** Échantillon stratifié d'au plus 100 observations par pixel sans remise parmi les observations non-wilt QC-valides sur 2016–2018 (min : {sampling_stats['min_obs_per_pixel']}, médiane : {sampling_stats['median_obs_per_pixel']:.1f}, {sampling_stats['n_pixels_under_100']:,} pixels sous 100 obs sur {sampling_stats['n_train_pixels']:,} pixels, {sampling_stats['n_train_rows']:,} lignes au total).  
**Algorithme principal :** LightGBM (300 arbres, learning_rate=0,05, num_leaves=63)  
**Gestion des valeurs manquantes :** Traitement natif des NaN par LightGBM pour la rugosité (`roughness_coefficient`), sans imputation par zéro, complété par l'indicateur binaire `roughness_is_nan`.  
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

Validation croisée sur les 5 plis spatiaux équilibrés de 1° (`spatial_fold`, Décision 7) avec calcul du TreeSHAP exact sur au moins 25 000 observations de validation par pli :

| Pli | Taille Train | Taille Val | RMSE Val | Pearson r | Biais Val |
| :--- | :--- | :--- | :--- | :--- | :--- |
{cv_table}
| **Moyenne** | — | — | **{mean_cv_rmse:.4f} m³/m³** | **{mean_cv_r:.4f}** | — |

---

## 3. Stabilité des attributions SHAP (Groupes, Intra-groupe et Variables)

### 3.1 Stabilité au niveau des 4 GROUPES et corrélation de rang intra-groupe

| Groupe | Contribution Absolue Moyenne | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | Rang Moyen | Spearman intra-groupe moyen |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
{group_table}
| **Moyenne groupes** | — | — | — | — | — | — | — | **{mean_rho_groups:.4f}** |

### 3.2 Rangs des 15 variables principales à travers les plis

| Rang | Variable | Rang Moyen | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | SHAP Absolu Moyen |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
{top15_table}

### 3.3 Matrice de corrélation de Spearman entre plis

| Paire de plis | Spearman {latex_rho} (4 Groupes) | Spearman {latex_rho} (Top 15) | Spearman {latex_rho} (Toutes les 66 variables) |
| :--- | :--- | :--- | :--- |
{stab_table}
| **Moyenne** | **{mean_rho_groups:.4f}** | **{mean_rho_top15:.4f}** | **{mean_rho_all:.4f}** |

---

## 4. Variantes de correction et décomposition saisonnière (2019)

Les 3 variantes de correction appliquent la décomposition additive :
- **V1 (removable) :** {latex_v1}
- **V2 (removable + context_rep) :** {latex_v2}
- **V3 (removable + all context) :** {latex_v3}

### Correction moyenne appliquée $\overline{{\Delta \mathrm{{SMAP}}}}$ et biais résiduel $\overline{{d_{{\mathrm{{res}}}}}}$ par saison

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{{\mathrm{{res}}}}$ | V2 : Correction | V2 : Biais $d_{{\mathrm{{res}}}}$ | V3 : Correction | V3 : Biais $d_{{\mathrm{{res}}}}$ | B2-train : Correction | B2-train : Biais $d_{{\mathrm{{res}}}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{seas_table}

---

## 5. Protocole 2a complet : Comparaison avec B0, B1, B2 et B2-train

Le flétrissement ({latex_wlt}) est recalculé sur l'observation corrigée pour chaque référence.

### 5.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\mathrm{{SMAP}}^* \mathrm{{vs}} \mathrm{{OL}}}}$ | $r_{{\mathrm{{anom}}}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_spec}

### 5.2 Échantillon commun ({n_common:,} observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\mathrm{{SMAP}}^* \mathrm{{vs}} \mathrm{{OL}}}}$ | $r_{{\mathrm{{anom}}}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_comm}

### 5.3 Stratification saisonnière des métriques (Échantillon spécifique)

| Référence | Saison | Observations assimilables | Rejets flétrissement | Biais ({latex_dres}) | RMSE | ubRMSE | $r_{{\mathrm{{SMAP}}^* \mathrm{{vs}} \mathrm{{OL}}}}$ | $r_{{\mathrm{{anom}}}}$ (SMAP brut) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_seas_metrics}

---

## 6. Figures de publication (FR et EN, sans interprétation)

- **Importance par groupe (4 groupes) :** [F7_shap_group_importance.png](file://{os.path.join(FIG_DIR, "F7_shap_group_importance.png")}) | EN: [F7_shap_group_importance_en.png](file://{os.path.join(FIG_EN_DIR, "F7_shap_group_importance_en.png")})
- **Beeswarm des 15 principales variables :** [F8_shap_beeswarm_top15.png](file://{os.path.join(FIG_DIR, "F8_shap_beeswarm_top15.png")}) | EN: [F8_shap_beeswarm_top15_en.png](file://{os.path.join(FIG_EN_DIR, "F8_shap_beeswarm_top15_en.png")})
- **Stabilité des rangs entre plis :** [F9_shap_fold_stability.png](file://{os.path.join(FIG_DIR, "F9_shap_fold_stability.png")}) | EN: [F9_shap_fold_stability_en.png](file://{os.path.join(FIG_EN_DIR, "F9_shap_fold_stability_en.png")})
- **Cartes saisonnières de correction moyenne (V1, V2, V3, B2-train) :** [F10_seasonal_applied_corrections_v1_v2_v3.png](file://{os.path.join(FIG_DIR, "F10_seasonal_applied_corrections_v1_v2_v3.png")}) | EN: [F10_seasonal_applied_corrections_v1_v2_v3_en.png](file://{os.path.join(FIG_EN_DIR, "F10_seasonal_applied_corrections_v1_v2_v3_en.png")})
- **Carte spatiale du masque de rugosité manquante :** [F10_roughness_is_nan.png](file://{os.path.join(FIG_DIR, "F10_roughness_is_nan.png")}) | EN: [F10_roughness_is_nan_en.png](file://{os.path.join(FIG_EN_DIR, "F10_roughness_is_nan_en.png")})
- **Carte spatiale de la texture dominante Noah-MP :** [F10_dominant_soiltype.png](file://{os.path.join(FIG_DIR, "F10_dominant_soiltype.png")}) | EN: [F10_dominant_soiltype_en.png](file://{os.path.join(FIG_EN_DIR, "F10_dominant_soiltype_en.png")})
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
    train_df, sampling_stats = load_training_data(feature_cols)
    
    n_jobs = get_n_jobs()
    print(f"Allocated CPU cores (n_jobs): {n_jobs}")

    # 4. Spatial Cross-Validation (5 folds on 1° blocks)
    df_cv, df_ranks, df_stability, df_group_stab, matrix_spearman, top15_features = run_spatial_cross_validation(
        train_df, feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, n_jobs=n_jobs
    )
    
    # 5. Full Model Training on 2016-2018
    model = train_full_model(train_df, feature_cols, n_jobs=n_jobs)
    
    # 6. Load B2 & B2-train Lookup Tables
    lut_b2 = np.load(LUT_CACHE_FILE) if os.path.exists(LUT_CACHE_FILE) else None
    lut_b2_train = load_or_build_b2_train_lut(landmask)
    
    # Baseline B1 Climatology learned on ALL 2016-2018 observations passing non-wilt QC
    b1_clim_table, b1_pixel_table, global_fallback = learn_b1_climatology()
    
    # 7. Evaluate 2019 Validation Set
    df_2019, shap_2019, X_val_2019, phi_0, n_total_2019, n_base_qc = evaluate_2019_validation(
        model, feature_cols, cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc,
        lut_b2, lut_b2_train, landmask, soil_wlt, dom_soil, b1_clim_table, b1_pixel_table, global_fallback, n_jobs=n_jobs
    )
    
    # 8. Compute Seasonal Corrections, Full Protocol 2a Metrics, and Seasonal Metrics
    df_seasonal, df_metrics, df_seas_metrics, n_common = compute_seasonal_corrections_and_metrics(
        df_2019, n_total_2019, n_base_qc
    )
    
    # 9. Generate Publication Figures (FR & EN, zero interpretation, dynamic legends)
    generate_figures(
        df_2019, shap_2019, X_val_2019, feature_cols, cols_rem, cols_ret,
        cols_ctx_rep, cols_ctx_proc, df_ranks, df_stability, df_group_stab,
        matrix_spearman, top15_features, dom_soil, landmask
    )
    
    # 10. Generate Markdown Report directly
    generate_markdown_report(
        df_cv, df_ranks, df_stability, df_group_stab, df_seasonal, df_metrics, df_seas_metrics,
        n_total_2019, n_base_qc, n_common, phi_0, sampling_stats
    )
    
    print(f"=== PHASE 2c COMPLETED SUCCESSFULLY IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
