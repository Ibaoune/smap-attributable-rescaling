#!/usr/bin/env python3
"""
Phase 4, Step 1: Evaluate 9 km EASE-grid reprojected V2 observations ("V2-9km") on 2019.
Reads back the generated 9 km SMAP HDF5 files from experiments/da_v2/obs/,
reprojects them to the NorthMor 5 km grid (replicating LIS neighbor_interp),
recomputes wilting point QC, and calculates all Décision 15 metrics for:
- B2_train_CDF
- LGB_V2_Rem_ContextRep (V2 at 5 km native modeling resolution)
- LGB_V2_9km (V2 at 9 km reprojected to 5 km)
on both Specific_Assimilable and Common_Assimilable populations.
Outputs comparison to reports/check_da_v2_9km_2019.csv.
No LightGBM retraining or model fitting.
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
import h5py
import netCDF4 as nc

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09, regrid_smap_to_northmor

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
CDF_B2_TRAIN_OBS = os.path.join(BASE_DIR, "experiments", "ldt_b2_train", "OUTPUT_obs", "smap_obs_cdf_NorthMor.nc")
LUT_B2_TRAIN_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_b2_train.npy")
CONTRIB_CACHE_FILE = os.path.join(BASE_DIR, "data", "interim", "lgb_contrib_2019.npy")
TABLE_2019 = os.path.join(BASE_DIR, "data", "processed", "modeling_table", "sar_modeling_table_2019.parquet")
OUT_OBS_DIR = os.path.join(BASE_DIR, "experiments", "da_v2", "obs")
CHECK_9KM_CSV = os.path.join(BASE_DIR, "reports", "check_da_v2_9km_2019.csv")

TEXTURE_NAMES = {
    1: 'Sand', 2: 'Loamy Sand', 3: 'Sandy Loam', 4: 'Silt Loam', 5: 'Silt',
    6: 'Loam', 7: 'Sandy Clay Loam', 8: 'Silty Clay Loam', 9: 'Clay Loam',
    10: 'Sandy Clay', 11: 'Silty Clay', 12: 'Clay', 13: 'Organic Material',
    14: 'Water', 15: 'Bedrock', 16: 'Other'
}


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
    t_start = time.time()
    print("=" * 80)
    print("EVALUATION DES OBSERVATIONS V2 REINTERPOLEES DEPUIS LE 9 KM (V2-9km)")
    print("ANNEE 2019 - COMPARAISON AVEC V2 (5 km) ET B2-train CDF")
    print("=" * 80)

    # 1. Grille et masques
    print(f"Chargement de la grille domaine 2D depuis {LIS_INPUT}...")
    with nc.Dataset(LIS_INPUT) as ds:
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
    row_map, col_map = latlon_to_ease2_m09(lat, lon)

    print(f"Chargement des masques et proprietes sol depuis {STATIC_FILE}...")
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]

    # 2. Chargement de la table de modelisation 2019
    print(f"Chargement de la table 2019 depuis {TABLE_2019}...")
    cols_rem, cols_ret, cols_ctx_rep, cols_ctx_proc, feature_cols = load_feature_groups()
    meta_cols = [
        'date', 'pass', 'pixel_id', 'month', 'season', 'smap_sm', 'ssm_ol_daily',
        'target_d', 'irrigfrac_gmia', 'gmia_class', 'spatial_block_id',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    load_cols = list(set(feature_cols + meta_cols))
    df_raw = pd.read_parquet(TABLE_2019, columns=load_cols)
    print(f"Lignes brutes 2019 : {len(df_raw):,}")

    qc_nw = (
        (df_raw['qc_lis_rain'] == 0) &
        (df_raw['qc_lis_frozen'] == 0) &
        (df_raw['qc_lis_forest'] == 0) &
        (df_raw['qc_lis_snow'] == 0) &
        (df_raw['qc_lis_high_veg'] == 0) &
        (df_raw['qc_lis_sat_limit'] == 0) &
        (df_raw['qc_lis_vwc_gt_5'] == 0) &
        (~df_raw['smap_sm'].isna())
    )
    df_2019 = df_raw[qc_nw].copy().reset_index(drop=True)
    n_base_qc = len(df_2019)
    print(f"Lignes valides QC non-wilt : {n_base_qc:,}")

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

    # 3. Modele B2-train CDF Matching
    print("Calcul de B2-train CDF matching...")
    lut_b2_train = np.load(LUT_B2_TRAIN_FILE)
    months = df_2019['month'].values
    with nc.Dataset(CDF_B2_TRAIN_OBS) as ds_o_tr:
        obs_x_tr = ds_o_tr.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta_tr = obs_x_tr[1, :, :] - obs_x_tr[0, :, :]
    obs_x0_tr = obs_x_tr[0, :, :]
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

    smap_star_b2_tr = apply_lut(lut_b2_train, obs_x0_tr, obs_delta_tr)
    df_2019['smap_star_b2_train'] = smap_star_b2_tr
    df_2019['d_res_b2_train'] = smap_star_b2_tr - df_2019['ssm_ol_daily']
    df_2019['under_wlt_b2_train'] = (df_2019['smap_star_b2_train'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b2_train'] = (
        (~df_2019['smap_star_b2_train'].isna()) &
        (df_2019['under_wlt_b2_train'] == 0) &
        (df_2019['smap_star_b2_train'] >= 0.01) &
        (df_2019['smap_star_b2_train'] <= 1.0)
    )

    # 4. Modele LGB_V2_Rem_ContextRep (V2 a 5 km)
    print("Calcul de LGB_V2_Rem_ContextRep (5 km) via le cache TreeSHAP...")
    contrib_2019 = np.load(CONTRIB_CACHE_FILE)
    phi_0 = float(contrib_2019[0, -1])
    shap_2019 = contrib_2019[:, :-1]
    idx_rem = [feature_cols.index(c) for c in cols_rem]
    idx_ctx_rep = [feature_cols.index(c) for c in cols_ctx_rep]
    phi_rem = shap_2019[:, idx_rem].sum(axis=1)
    phi_ctx_rep = shap_2019[:, idx_ctx_rep].sum(axis=1)

    df_2019['smap_star_v2'] = df_2019['smap_sm'] - (phi_rem + phi_ctx_rep + phi_0)
    df_2019['d_res_v2'] = df_2019['smap_star_v2'] - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v2'] = (df_2019['smap_star_v2'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v2'] = (
        (df_2019['under_wlt_v2'] == 0) &
        (df_2019['smap_star_v2'] >= 0.01) &
        (df_2019['smap_star_v2'] <= 1.0) &
        (~df_2019['smap_star_v2'].isna())
    )

    # 5. Relecture des fichiers 9 km et reinterpolation a 5 km (V2-9km)
    print("Relecture des fichiers 9 km HDF5 et projection 5 km (comme LIS)...")
    raw_files = sorted(glob.glob(os.path.join(OUT_OBS_DIR, "*", "*.h5")))
    print(f"Fichiers HDF5 trouves dans {OUT_OBS_DIR} : {len(raw_files)}")

    df_2019['row_in_df'] = np.arange(len(df_2019))
    groups_for_read = {k: v for k, v in df_2019.groupby(['date', 'pass'])}
    smap_star_v2_9km = np.full(n_base_qc, np.nan, dtype=np.float32)

    for i, h5_path in enumerate(raw_files):
        fname = os.path.basename(h5_path)
        date_str = fname.split('_')[5]
        yyyy, mm, dd = date_str[:4], date_str[4:6], date_str[6:8]
        date_iso = f"{yyyy}-{mm}-{dd}"

        with h5py.File(h5_path, 'r') as h5_in:
            for p_code, grp_name, sm_var in [
                ('AM', 'Soil_Moisture_Retrieval_Data_AM', 'soil_moisture'),
                ('PM', 'Soil_Moisture_Retrieval_Data_PM', 'soil_moisture_pm')
            ]:
                if grp_name not in h5_in:
                    continue

                key_gp = (date_iso, p_code)
                if key_gp in groups_for_read:
                    sub = groups_for_read[key_gp]
                    sm_arr_read = h5_in[grp_name][sm_var][:]
                    regridded_domain = regrid_smap_to_northmor(sm_arr_read, row_map, col_map).ravel()
                    vals = regridded_domain[sub['pixel_id'].values]
                    # Filter out missing values
                    vals_clean = np.where(vals <= -9999.0, np.nan, vals)
                    smap_star_v2_9km[sub['row_in_df'].values] = vals_clean

        if (i + 1) % 50 == 0 or (i + 1) == len(raw_files):
            print(f"  [{i+1}/{len(raw_files)}] Fichiers relus et reinterpoles...")

    df_2019['smap_star_v2_9km'] = smap_star_v2_9km
    df_2019['d_res_v2_9km'] = smap_star_v2_9km - df_2019['ssm_ol_daily']
    df_2019['under_wlt_v2_9km'] = (df_2019['smap_star_v2_9km'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_v2_9km'] = (
        (~df_2019['smap_star_v2_9km'].isna()) &
        (df_2019['under_wlt_v2_9km'] == 0) &
        (df_2019['smap_star_v2_9km'] >= 0.01) &
        (df_2019['smap_star_v2_9km'] <= 1.0)
    )

    # 6. Echantillon commun (Common_Assimilable)
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b2_train'] &
        df_2019['assimilable_v2'] &
        df_2019['assimilable_v2_9km']
    )
    print(f"Echantillon commun assimilable (B2-train & V2 & V2-9km) : {df_2019['assimilable_common'].sum():,} observations")

    # 7. Calcul de toutes les metriques de la Decision 15
    models = [
        ('B2_train_CDF', 'd_res_b2_train', 'smap_star_b2_train', 'assimilable_b2_train', 'under_wlt_b2_train'),
        ('LGB_V2_Rem_ContextRep', 'd_res_v2', 'smap_star_v2', 'assimilable_v2', 'under_wlt_v2'),
        ('LGB_V2_9km', 'd_res_v2_9km', 'smap_star_v2_9km', 'assimilable_v2_9km', 'under_wlt_v2_9km')
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

    df_out = pd.DataFrame(metrics_rows)
    df_out.to_csv(CHECK_9KM_CSV, index=False)
    print(f"\nResultats sauvegardes dans {CHECK_9KM_CSV}")

    print("\n" + "=" * 100)
    print("TABLEAU COMPARATIF DECISION 15 : B2-train vs V2 (5 km) vs V2-9km (2019)")
    print("=" * 100)
    cols_display = [
        'evaluation_population', 'model', 'n_assimilable', 'assimilable_pct',
        'wilt_rejection_pct', 'bias', 'rmse', 'r_smapstar_vs_ol', 'r_anom_with_smap_raw',
        'gmia_contrast_cl_jja'
    ]
    print(df_out[cols_display].to_string(index=False))
    print(f"\nExecution terminee en {time.time()-t_start:.1f}s.")


if __name__ == '__main__':
    main()
