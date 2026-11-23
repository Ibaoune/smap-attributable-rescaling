#!/usr/bin/env python3
"""
scripts/07_run_phase2a_baselines.py

PHASE 2a : Baseline benchmark before machine learning models.
Strict protocol:
- Training: 2016-2018 (split == 'train')
- Evaluation: 2019 (split == 'val')
- 2020 remains strictly sealed!

Protocol refinements (Point de contrôle 5b):
1. Evaluation population: start from observations passing all QC rules EXCEPT wilt limit.
   Apply correction, then recalculate the wilt limit rule on corrected SMAP (smap* < theta_wlt + 0.02).
   Metrics (Bias, RMSE) computed on assimilable observations specific to each reference,
   and on the common assimilable sample across all three references.
2. B1 training: learned on 2016-2018 train split on all observations passing non-wilt QC rules.
3. B2 LDT metadata: documented period (2016-2021) in report header, unrescaled counts reported.
4. Additional metrics: assimilable shares, GMIA >30% vs 0% contrast with 1° block bootstrap CIs,
   anomaly correlation with raw SMAP.
5. Relative BASE_DIR.
6. Zero interpretive phrasing (only computed facts).
"""

import os
import glob
import time
import netCDF4 as nc
import numpy as np
import pandas as pd

# Paths (Relative to repository root)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
CDF_MODEL_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "noahmp_model_cdf_NorthMor.nc")
CDF_OBS_FILE = os.path.join(BASE_DIR, "data", "raw_copied", "smap_obs_cdf_NorthMor.nc")
LUT_CACHE_FILE = os.path.join(BASE_DIR, "data", "interim", "lut_cdf_rescaling.npy")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
FIG_DIR = os.path.join(BASE_DIR, "figures", "descriptive")

DA_CDF_BASE = os.path.join(
    os.path.dirname(BASE_DIR),
    "NoahMP_Morocco",
    "experiments",
    "NorthMor",
    "matrix_2016_2020",
    "DA_cdf_noirr_2016_2020",
    "output"
)

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def load_or_build_cdf_lut(landmask):
    """Load cached CDF rescaling LUT or build it from raw netCDF files."""
    if os.path.exists(LUT_CACHE_FILE):
        print(f"Loading cached CDF LUT from {LUT_CACHE_FILE}...")
        lut = np.load(LUT_CACHE_FILE)
        print(f"Loaded LUT: shape {lut.shape}, valid entries: {(~np.isnan(lut)).sum():,}")
        return lut

    print("Building CDF rescaling lookup table from LDT files...")
    t0 = time.time()
    with nc.Dataset(CDF_MODEL_FILE) as ds_m:
        model_x = ds_m.variables['SoilMoist_xrange'][:, 0, :, :] # (100, 12, 16030)
        model_cdf = ds_m.variables['SoilMoist_CDF'][:, 0, :, :]
    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
        obs_cdf = ds_o.variables['SoilMoist_CDF'][:, 0, :, :]

    nbins = 100
    ntimes = 12
    ngp = 16030
    lut = np.full((nbins, ntimes, ngp), np.nan, dtype=np.float32)

    for m in range(ntimes):
        t_m0 = time.time()
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
        print(f"  Month {m+1:02d}/12 processed in {time.time()-t_m0:.1f}s")

    os.makedirs(os.path.dirname(LUT_CACHE_FILE), exist_ok=True)
    np.save(LUT_CACHE_FILE, lut)
    print(f"CDF LUT built and saved in {time.time()-t0:.1f}s: {LUT_CACHE_FILE}")
    return lut


def verify_cdf_matching_against_daobs(lut, landmask):
    """Verify offline CDF matching against DAOBS Record 2 for 5 sample dates."""
    print("--- Verifying CDF Matching against DAOBS Record 2 ---")
    tile_map = np.full((121, 181), -1, dtype=int)
    t = 0
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                tile_map[r, c] = t
                t += 1

    sample_specs = [
        ("2016-01-01 06:00", os.path.join(DA_CDF_BASE, "2016-01/DAOBS/201601/LISDAOBS_201601010600.a01.d01.1gs4r"), 0),
        ("2019-01-01 06:00", os.path.join(DA_CDF_BASE, "2019-01/DAOBS/201901/LISDAOBS_201901010600.a01.d01.1gs4r"), 0),
        ("2019-04-01 06:00", os.path.join(DA_CDF_BASE, "2019-04/DAOBS/201904/LISDAOBS_201904010600.a01.d01.1gs4r"), 3),
        ("2019-07-24 06:00", os.path.join(DA_CDF_BASE, "2019-07/DAOBS/201907/LISDAOBS_201907240600.a01.d01.1gs4r"), 6),
        ("2019-10-01 18:00", os.path.join(DA_CDF_BASE, "2019-10/DAOBS/201910/LISDAOBS_201910011800.a01.d01.1gs4r"), 9),
    ]

    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta = obs_x[1, :, :] - obs_x[0, :, :] # (12, 16030)
    obs_x0 = obs_x[0, :, :]

    verification_results = []
    for label, fpath, m_idx in sample_specs:
        if not os.path.exists(fpath):
            print(f"  Warning: DAOBS file not found: {fpath}")
            continue

        with open(fpath, "rb") as f:
            _ = np.fromfile(f, dtype=">i4", count=1)
            rec1 = np.fromfile(f, dtype=">f4", count=21901).reshape(121, 181)
            _ = np.fromfile(f, dtype=">i4", count=2)
            rec2 = np.fromfile(f, dtype=">f4", count=21901).reshape(121, 181)

        diffs = []
        for r in range(121):
            for c in range(181):
                tile = tile_map[r, c]
                if tile == -1:
                    continue
                val = rec1[r, c]
                if val <= -9999.0 or np.isnan(val) or val < 0.01 or val > 1.0:
                    continue
                daobs_val = rec2[r, c]
                if daobs_val <= 0.01 or daobs_val > 1.0:
                    continue

                d_val = obs_delta[m_idx, tile]
                if d_val <= 0:
                    continue
                binval = int(np.floor((val - obs_x0[m_idx, tile]) / d_val + 0.5)) + 1
                binval = max(1, min(100, binval)) - 1
                m_val = lut[binval, m_idx, tile]
                if not np.isnan(m_val):
                    diffs.append(abs(m_val - daobs_val))

        diffs = np.array(diffs)
        n_pts = len(diffs)
        max_diff = float(diffs.max()) if n_pts > 0 else np.nan
        mean_diff = float(diffs.mean()) if n_pts > 0 else np.nan
        pct_exact = float((diffs < 1e-5).mean() * 100) if n_pts > 0 else np.nan
        print(f"  {label}: {n_pts:,} valid points, max abs diff: {max_diff:.8f}, mean abs diff: {mean_diff:.8f}, exact: {pct_exact:.2f}%")
        verification_results.append({
            "sample_date": label,
            "n_points": n_pts,
            "mean_abs_diff": mean_diff,
            "max_abs_diff": max_diff,
            "pct_exact_match": pct_exact
        })

    return verification_results


def learn_b1_climatology():
    """Learn B1 pixel-month mean bias on 2016-2018 (split == 'train') on observations passing non-wilt QC."""
    print("--- Learning B1 Climatology on 2016-2018 (Train Split, Non-wilt QC) ---")
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
    print(f"Loaded {len(train_df):,} training observations (non-wilt QC valid, 2016-2018)")

    climatology = train_df.groupby(['pixel_id', 'month'])['target_d'].mean().to_dict()
    pixel_fallback = train_df.groupby('pixel_id')['target_d'].mean().to_dict()
    global_fallback = float(train_df['target_d'].mean())
    print(f"Climatology learned: {len(climatology):,} (pixel, month) cells, global mean d: {global_fallback:.6f}")
    return climatology, pixel_fallback, global_fallback


def evaluate_baselines_on_2019(lut, climatology, pixel_fallback, global_fallback, landmask, soil_wlt, dom_soil):
    """Evaluate B0, B1, B2 on 2019 validation set starting from non-wilt QC rules."""
    print("--- Evaluating Baselines on 2019 (Val Split) ---")
    pq_2019 = os.path.join(TABLE_DIR, "sar_modeling_table_2019.parquet")
    cols = [
        'pixel_id', 'month', 'season', 'gmia_class', 'spatial_block_id',
        'smap_sm', 'ssm_ol_daily', 'target_d', 'split',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    df_raw_2019 = pd.read_parquet(pq_2019, columns=cols)
    n_total_2019 = len(df_raw_2019)
    print(f"Loaded 2019 table: {n_total_2019:,} total rows")

    # Filter base population: all observations passing non-wilt QC rules
    qc_nw = (
        (df_raw_2019['qc_lis_rain'] == 0) &
        (df_raw_2019['qc_lis_frozen'] == 0) &
        (df_raw_2019['qc_lis_forest'] == 0) &
        (df_raw_2019['qc_lis_snow'] == 0) &
        (df_raw_2019['qc_lis_high_veg'] == 0) &
        (df_raw_2019['qc_lis_sat_limit'] == 0) &
        (df_raw_2019['qc_lis_vwc_gt_5'] == 0)
    )
    df_2019 = df_raw_2019[qc_nw].copy()
    n_base_qc = len(df_2019)
    print(f"Base population (non-wilt QC valid): {n_base_qc:,} rows ({n_base_qc/n_total_2019*100:.2f}%)")

    # Tile map
    tile_map = np.full((121, 181), -1, dtype=int)
    t = 0
    for r in range(121):
        for c in range(181):
            if landmask[r, c] == 1:
                tile_map[r, c] = t
                t += 1

    # Extract pixel metadata arrays
    r_idx = df_2019['pixel_id'].values // 181
    c_idx = df_2019['pixel_id'].values % 181
    tiles = tile_map[r_idx, c_idx]
    wlt_thresh = soil_wlt[r_idx, c_idx] + 0.02
    soil_ids = dom_soil[r_idx, c_idx]
    soil_names = [TEXTURE_NAMES.get(s, "Other") for s in soil_ids]
    df_2019['texture_name'] = soil_names
    df_2019['wlt_thresh'] = wlt_thresh

    # Baseline B0: No correction
    df_2019['smap_star_b0'] = df_2019['smap_sm']
    df_2019['d_res_b0'] = df_2019['target_d']
    df_2019['under_wlt_b0'] = (df_2019['smap_star_b0'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b0'] = (
        (df_2019['under_wlt_b0'] == 0) &
        (df_2019['smap_star_b0'] >= 0.01) &
        (df_2019['smap_star_b0'] <= 1.0) &
        (~df_2019['smap_star_b0'].isna())
    )

    # Baseline B1: Climatology
    p_ids = df_2019['pixel_id'].values
    months = df_2019['month'].values
    c_hat = np.array([climatology.get((p, m), pixel_fallback.get(p, global_fallback)) for p, m in zip(p_ids, months)], dtype=np.float32)
    df_2019['c_hat_b1'] = c_hat
    df_2019['smap_star_b1'] = df_2019['smap_sm'] - c_hat
    df_2019['d_res_b1'] = df_2019['target_d'] - c_hat
    df_2019['under_wlt_b1'] = (df_2019['smap_star_b1'] < df_2019['wlt_thresh']).astype(int)
    df_2019['assimilable_b1'] = (
        (df_2019['under_wlt_b1'] == 0) &
        (df_2019['smap_star_b1'] >= 0.01) &
        (df_2019['smap_star_b1'] <= 1.0) &
        (~df_2019['smap_star_b1'].isna())
    )

    # Baseline B2: CDF v2 offline
    with nc.Dataset(CDF_OBS_FILE) as ds_o:
        obs_x = ds_o.variables['SoilMoist_xrange'][:, 0, :, :]
    obs_delta = obs_x[1, :, :] - obs_x[0, :, :] # (12, 16030)
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

    # Common sample (assimilable in B0, B1, and B2 simultaneously)
    df_2019['assimilable_common'] = (
        df_2019['assimilable_b0'] &
        df_2019['assimilable_b1'] &
        df_2019['assimilable_b2']
    )

    return df_2019, n_total_2019, n_base_qc


def compute_b2_missing_stats(df_2019):
    """Compute count and percentage where B2 returns no value (IQR failure, out of bounds), overall and by season."""
    print("--- Computing B2 Missing/Unrescaled Statistics ---")
    rows = []
    seasons = ["All", "DJF", "MAM", "JJA", "SON"]
    for s in seasons:
        sub = df_2019 if s == "All" else df_2019[df_2019['season'] == s]
        n_sub = len(sub)
        n_missing = int((sub['b2_has_value'] == 0).sum())
        pct_missing = (n_missing / n_sub * 100.0) if n_sub > 0 else 0.0
        n_valid = n_sub - n_missing
        pct_valid = (n_valid / n_sub * 100.0) if n_sub > 0 else 0.0
        rows.append({
            "season": s,
            "total_obs_no_wilt": n_sub,
            "b2_valid_count": n_valid,
            "b2_valid_pct": pct_valid,
            "b2_missing_count": n_missing,
            "b2_missing_pct": pct_missing
        })
    return pd.DataFrame(rows)


def compute_anomaly_correlations(df_2019):
    """Compute anomaly correlation between corrected SMAP and raw SMAP on specific and common assimilable samples."""
    print("--- Computing Anomaly Correlations ---")
    results = {}

    for sample_name, mask in [('specific', None), ('common', df_2019['assimilable_common'])]:
        sub_b0 = df_2019[df_2019['assimilable_b0']] if mask is None else df_2019[mask]
        sub_b1 = df_2019[df_2019['assimilable_b1']] if mask is None else df_2019[mask]
        sub_b2 = df_2019[df_2019['assimilable_b2']] if mask is None else df_2019[mask]

        # B0: raw SMAP vs raw SMAP
        anom_raw_b0 = sub_b0['smap_sm'] - sub_b0.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')
        anom_star_b0 = sub_b0['smap_star_b0'] - sub_b0.groupby(['pixel_id', 'month'])['smap_star_b0'].transform('mean')
        r_b0 = float(np.corrcoef(anom_raw_b0, anom_star_b0)[0, 1])

        # B1: corrected B1 vs raw SMAP
        anom_raw_b1 = sub_b1['smap_sm'] - sub_b1.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')
        anom_star_b1 = sub_b1['smap_star_b1'] - sub_b1.groupby(['pixel_id', 'month'])['smap_star_b1'].transform('mean')
        r_b1 = float(np.corrcoef(anom_raw_b1, anom_star_b1)[0, 1])

        # B2: corrected B2 vs raw SMAP
        anom_raw_b2 = sub_b2['smap_sm'] - sub_b2.groupby(['pixel_id', 'month'])['smap_sm'].transform('mean')
        anom_star_b2 = sub_b2['smap_star_b2'] - sub_b2.groupby(['pixel_id', 'month'])['smap_star_b2'].transform('mean')
        r_b2 = float(np.corrcoef(anom_raw_b2, anom_star_b2)[0, 1])

        results[sample_name] = {
            'r_anom_b0': r_b0,
            'r_anom_b1': r_b1,
            'r_anom_b2': r_b2,
            'n_b0': len(sub_b0),
            'n_b1': len(sub_b1),
            'n_b2': len(sub_b2)
        }
    return results


def compute_gmia_contrast_with_bootstrap(df_2019):
    """Compute GMIA >30% vs GMIA 0% contrast for residual d with 1° block bootstrap 95% CIs and distinct pixel counts."""
    print("--- Computing GMIA Contrast (>30% vs 0%) at Fixed Texture and Season ---")
    np.random.seed(42)
    B = 500
    rows = []

    # Focus on Clay Loam (balanced irrigation representation) and Loam (dominant regional texture)
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

            for ref, col_res, mask_col in [
                ('b0', 'd_res_b0', 'assimilable_b0' if sample_type == 'specific' else 'assimilable_common'),
                ('b1', 'd_res_b1', 'assimilable_b1' if sample_type == 'specific' else 'assimilable_common'),
                ('b2', 'd_res_b2', 'assimilable_b2' if sample_type == 'specific' else 'assimilable_common')
            ]:
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

                # Block bootstrap on 1° spatial blocks
                sub_both = sub_ref[sub_ref['gmia_class'].isin(['0%', '>30%'])]
                n_blocks_total = int(sub_both['spatial_block_id'].nunique())
                if n_obs_0 >= 10 and n_obs_gt30 >= 10 and sub_both['spatial_block_id'].nunique() > 1:
                    stats = sub_both.groupby(['spatial_block_id', 'gmia_class'])[col_res].agg(['sum', 'count']).unstack('gmia_class').fillna(0)
                    sum_0_arr = stats[('sum', '0%')].values if ('sum', '0%') in stats else np.zeros(len(stats))
                    cnt_0_arr = stats[('count', '0%')].values if ('count', '0%') in stats else np.zeros(len(stats))
                    sum_gt30_arr = stats[('sum', '>30%')].values if ('sum', '>30%') in stats else np.zeros(len(stats))
                    cnt_gt30_arr = stats[('count', '>30%')].values if ('count', '>30%') in stats else np.zeros(len(stats))
                    n_blocks = len(stats)

                    boot_diffs = np.empty(B)
                    for b_i in range(B):
                        idx = np.random.choice(n_blocks, size=n_blocks, replace=True)
                        c0 = cnt_0_arr[idx].sum()
                        cgt = cnt_gt30_arr[idx].sum()
                        m0 = (sum_0_arr[idx].sum() / c0) if c0 > 0 else np.nan
                        mgt = (sum_gt30_arr[idx].sum() / cgt) if cgt > 0 else np.nan
                        boot_diffs[b_i] = mgt - m0
                    valid_boot = boot_diffs[~np.isnan(boot_diffs)]
                    if len(valid_boot) > 50:
                        ci_low = float(np.percentile(valid_boot, 2.5))
                        ci_high = float(np.percentile(valid_boot, 97.5))
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

    df_contrast = pd.DataFrame(rows)
    return df_contrast


def compute_comprehensive_strata_metrics(df_2019, n_total_2019, n_base_qc):
    """Compute Bias and RMSE on both specific assimilable sets and common sample across strata."""
    print("--- Computing Comprehensive Stratified Metrics (Specific & Common) ---")
    results = []

    top_textures = ["Loam", "Clay Loam", "Sandy Clay Loam", "Sandy Loam", "Clay"]
    seasons = ["DJF", "MAM", "JJA", "SON"]
    gmia_classes = ["0%", "0-10%", "10-30%", ">30%"]

    def eval_group(sub_df, stratum_type, stratum_name):
        n_sub_base = len(sub_df)

        # 1. Specific assimilable sets
        s_b0 = sub_df[sub_df['assimilable_b0']]
        s_b1 = sub_df[sub_df['assimilable_b1']]
        s_b2 = sub_df[sub_df['assimilable_b2']]

        n_b0 = len(s_b0)
        n_b1 = len(s_b1)
        n_b2 = len(s_b2)

        bias_spec_b0 = float(s_b0['d_res_b0'].mean()) if n_b0 > 0 else np.nan
        rmse_spec_b0 = float(np.sqrt((s_b0['d_res_b0']**2).mean())) if n_b0 > 0 else np.nan

        bias_spec_b1 = float(s_b1['d_res_b1'].mean()) if n_b1 > 0 else np.nan
        rmse_spec_b1 = float(np.sqrt((s_b1['d_res_b1']**2).mean())) if n_b1 > 0 else np.nan

        bias_spec_b2 = float(s_b2['d_res_b2'].mean()) if n_b2 > 0 else np.nan
        rmse_spec_b2 = float(np.sqrt((s_b2['d_res_b2']**2).mean())) if n_b2 > 0 else np.nan

        # 2. Common assimilable sample
        s_com = sub_df[sub_df['assimilable_common']]
        n_com = len(s_com)

        bias_com_b0 = float(s_com['d_res_b0'].mean()) if n_com > 0 else np.nan
        rmse_com_b0 = float(np.sqrt((s_com['d_res_b0']**2).mean())) if n_com > 0 else np.nan

        bias_com_b1 = float(s_com['d_res_b1'].mean()) if n_com > 0 else np.nan
        rmse_com_b1 = float(np.sqrt((s_com['d_res_b1']**2).mean())) if n_com > 0 else np.nan

        bias_com_b2 = float(s_com['d_res_b2'].mean()) if n_com > 0 else np.nan
        rmse_com_b2 = float(np.sqrt((s_com['d_res_b2']**2).mean())) if n_com > 0 else np.nan

        return {
            'stratum_type': stratum_type,
            'stratum_name': stratum_name,
            'n_base_qc': n_sub_base,
            # Specific
            'n_b0_spec': n_b0,
            'pct_b0_spec': (n_b0 / n_sub_base * 100.0) if n_sub_base > 0 else np.nan,
            'bias_b0_spec': bias_spec_b0,
            'rmse_b0_spec': rmse_spec_b0,
            'n_b1_spec': n_b1,
            'pct_b1_spec': (n_b1 / n_sub_base * 100.0) if n_sub_base > 0 else np.nan,
            'bias_b1_spec': bias_spec_b1,
            'rmse_b1_spec': rmse_spec_b1,
            'n_b2_spec': n_b2,
            'pct_b2_spec': (n_b2 / n_sub_base * 100.0) if n_sub_base > 0 else np.nan,
            'bias_b2_spec': bias_spec_b2,
            'rmse_b2_spec': rmse_spec_b2,
            # Common
            'n_common': n_com,
            'pct_common': (n_com / n_sub_base * 100.0) if n_sub_base > 0 else np.nan,
            'bias_b0_com': bias_com_b0,
            'rmse_b0_com': rmse_com_b0,
            'bias_b1_com': bias_com_b1,
            'rmse_b1_com': rmse_com_b1,
            'bias_b2_com': bias_com_b2,
            'rmse_b2_com': rmse_com_b2,
        }

    # Global
    results.append(eval_group(df_2019, "Overall", "Global (All)"))

    # Season
    for s in seasons:
        results.append(eval_group(df_2019[df_2019['season'] == s], "Season", s))

    # Dominant Texture (Top 5)
    for tex in top_textures:
        results.append(eval_group(df_2019[df_2019['texture_name'] == tex], "Texture", tex))

    # GMIA Class
    for gmia in gmia_classes:
        results.append(eval_group(df_2019[df_2019['gmia_class'] == gmia], "GMIA", gmia))

    df_metrics = pd.DataFrame(results)
    return df_metrics


def generate_refined_synthesis_report(metrics_df, b2_missing_df, contrast_df, anom_corr, daobs_verif, n_total_2019, n_base_qc):
    """Generate Markdown report dynamically from computed metrics with ZERO interpretive words."""
    print("--- Generating Refined Phase 2a Synthesis Report (Zero Interpretive Phrasing) ---")
    out_csv = os.path.join(REPORTS_DIR, "phase2a_baselines_metrics.csv")
    metrics_df.to_csv(out_csv, index=False)
    print(f"Saved metrics CSV: {out_csv}")

    contrast_csv = os.path.join(REPORTS_DIR, "phase2a_gmia_contrast.csv")
    contrast_df.to_csv(contrast_csv, index=False)
    print(f"Saved contrast CSV: {contrast_csv}")

    glob = metrics_df[metrics_df['stratum_type'] == "Overall"].iloc[0]

    # DAOBS table rows
    daobs_rows = []
    for d in daobs_verif:
        daobs_rows.append(
            f"| `{d['sample_date']}` | {d['n_points']:,} | {d['mean_abs_diff']:.8f} | {d['max_abs_diff']:.8f} | {d['pct_exact_match']:.2f}% |"
        )
    daobs_table = "\n".join(daobs_rows)

    # B2 Missing Table rows
    b2_miss_rows = []
    for _, row in b2_missing_df.iterrows():
        b2_miss_rows.append(
            f"| {row['season']} | {row['total_obs_no_wilt']:,} | {row['b2_valid_count']:,} | {row['b2_valid_pct']:.2f}% | {row['b2_missing_count']:,} | {row['b2_missing_pct']:.2f}% |"
        )
    b2_miss_table = "\n".join(b2_miss_rows)

    # Assimilable population table
    pop_table = f"""| Échantillon | Référence | Observations assimilables | Part du socle hors flétrissement ({n_base_qc:,} obs) | Part du total 2019 ({n_total_2019:,} obs) |
|:---|:---|:---:|:---:|:---:|
| Spécifique | **B0 (NoCDF)** | {glob['n_b0_spec']:,} | {glob['pct_b0_spec']:.2f}% | {glob['n_b0_spec']/n_total_2019*100.0:.2f}% |
| Spécifique | **B1 (Climatologie)** | {glob['n_b1_spec']:,} | {glob['pct_b1_spec']:.2f}% | {glob['n_b1_spec']/n_total_2019*100.0:.2f}% |
| Spécifique | **B2 (CDF v2)** | {glob['n_b2_spec']:,} | {glob['pct_b2_spec']:.2f}% | {glob['n_b2_spec']/n_total_2019*100.0:.2f}% |
| Commun | **B0 / B1 / B2** | {glob['n_common']:,} | {glob['pct_common']:.2f}% | {glob['n_common']/n_total_2019*100.0:.2f}% |"""

    # Global comparative metrics table
    glob_metrics_table = f"""| Référence | Biais résiduel (Échantillon spécifique) | RMSE résiduelle (Échantillon spécifique) | Biais résiduel (Échantillon commun) | RMSE résiduelle (Échantillon commun) |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | {glob['bias_b0_spec']:+.4f} m³/m³ | {glob['rmse_b0_spec']:.4f} m³/m³ | {glob['bias_b0_com']:+.4f} m³/m³ | {glob['rmse_b0_com']:.4f} m³/m³ |
| **B1 (Climatologie)** | {glob['bias_b1_spec']:+.4f} m³/m³ | {glob['rmse_b1_spec']:.4f} m³/m³ | {glob['bias_b1_com']:+.4f} m³/m³ | {glob['rmse_b1_com']:.4f} m³/m³ |
| **B2 (CDF v2)** | {glob['bias_b2_spec']:+.4f} m³/m³ | {glob['rmse_b2_spec']:.4f} m³/m³ | {glob['bias_b2_com']:+.4f} m³/m³ | {glob['rmse_b2_com']:.4f} m³/m³ |"""

    # Season table (specific & common)
    season_rows = []
    for s in ["DJF", "MAM", "JJA", "SON"]:
        s_row = metrics_df[metrics_df['stratum_name'] == s].iloc[0]
        season_rows.append(
            f"| {s} | {s_row['rmse_b0_spec']:.4f} | {s_row['rmse_b1_spec']:.4f} | {s_row['rmse_b2_spec']:.4f} | "
            f"{s_row['n_common']:,} | {s_row['rmse_b0_com']:.4f} | {s_row['rmse_b1_com']:.4f} | {s_row['rmse_b2_com']:.4f} |"
        )
    season_table = "\n".join(season_rows)

    # Texture table (specific & common)
    tex_rows = []
    for tex in ["Loam", "Clay Loam", "Sandy Clay Loam", "Sandy Loam", "Clay"]:
        t_row = metrics_df[metrics_df['stratum_name'] == tex].iloc[0]
        tex_rows.append(
            f"| {tex} | {t_row['n_b0_spec']:,} | {t_row['rmse_b0_spec']:.4f} | "
            f"{t_row['n_b1_spec']:,} | {t_row['rmse_b1_spec']:.4f} | "
            f"{t_row['n_b2_spec']:,} | {t_row['rmse_b2_spec']:.4f} | "
            f"{t_row['n_common']:,} | {t_row['rmse_b1_com']:.4f} | {t_row['rmse_b2_com']:.4f} |"
        )
    tex_table = "\n".join(tex_rows)

    # GMIA table (specific & common)
    gmia_rows = []
    for g in ["0%", "0-10%", "10-30%", ">30%"]:
        g_row = metrics_df[metrics_df['stratum_name'] == g].iloc[0]
        gmia_rows.append(
            f"| {g} | {g_row['n_b0_spec']:,} | {g_row['rmse_b0_spec']:.4f} | "
            f"{g_row['n_b1_spec']:,} | {g_row['rmse_b1_spec']:.4f} | "
            f"{g_row['n_b2_spec']:,} | {g_row['rmse_b2_spec']:.4f} | "
            f"{g_row['n_common']:,} | {g_row['rmse_b1_com']:.4f} | {g_row['rmse_b2_com']:.4f} |"
        )
    gmia_table = "\n".join(gmia_rows)

    # Anomaly correlation table
    anom_spec = anom_corr['specific']
    anom_com = anom_corr['common']
    anom_table = f"""| Référence | Anomaly Correlation (Échantillon spécifique) | N observations spécifique | Anomaly Correlation (Échantillon commun) | N observations commun |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | {anom_spec['r_anom_b0']:.4f} | {anom_spec['n_b0']:,} | {anom_com['r_anom_b0']:.4f} | {anom_com['n_b0']:,} |
| **B1 (Climatologie)** | {anom_spec['r_anom_b1']:.4f} | {anom_spec['n_b1']:,} | {anom_com['r_anom_b1']:.4f} | {anom_com['n_b1']:,} |
| **B2 (CDF v2)** | {anom_spec['r_anom_b2']:.4f} | {anom_spec['n_b2']:,} | {anom_com['r_anom_b2']:.4f} | {anom_com['n_b2']:,} |"""

    # GMIA Contrast Table (Clay Loam across seasons, specific and common)
    contrast_rows = []
    for _, c_row in contrast_df.iterrows():
        contrast_rows.append(
            f"| {c_row['texture']} | {c_row['season']} | {c_row['sample_type']} | "
            f"{c_row['b0_contrast']:+.4f} [{c_row['b0_ci_low']:+.4f}, {c_row['b0_ci_high']:+.4f}] ({c_row['b0_n_px_gt30']}/{c_row['b0_n_px_0']} px, {c_row['b0_n_blocks_gt30']}/{c_row['b0_n_blocks_0']} blk) | "
            f"{c_row['b1_contrast']:+.4f} [{c_row['b1_ci_low']:+.4f}, {c_row['b1_ci_high']:+.4f}] ({c_row['b1_n_px_gt30']}/{c_row['b1_n_px_0']} px, {c_row['b1_n_blocks_gt30']}/{c_row['b1_n_blocks_0']} blk) | "
            f"{c_row['b2_contrast']:+.4f} [{c_row['b2_ci_low']:+.4f}, {c_row['b2_ci_high']:+.4f}] ({c_row['b2_n_px_gt30']}/{c_row['b2_n_px_0']} px, {c_row['b2_n_blocks_gt30']}/{c_row['b2_n_blocks_0']} blk) |"
        )
    contrast_table = "\n".join(contrast_rows)

    # Load F3 table for inclusion
    f3_csv_path = os.path.join(FIG_DIR, "F3_d_by_gmia_texture_season.csv")
    f3_table_rows = []
    if os.path.exists(f3_csv_path):
        f3_df = pd.read_csv(f3_csv_path)
        for _, r_f3 in f3_df.iterrows():
            f3_table_rows.append(
                f"| {r_f3['soil_name']} | {r_f3['season_group_fr']} | {r_f3['gmia_class']} | "
                f"{int(r_f3['n_obs']):,} | {int(r_f3['n_distinct_pixels']):,} | {float(r_f3['mean']):+.4f} ± {float(r_f3['ci95_half']):.4f} |"
            )
    f3_table = "\n".join(f3_table_rows[:25]) # Top 25 rows in main text

    # Monthly wilt diagnostic from F5
    f5_csv_path = os.path.join(FIG_DIR, "F5_wilt_diagnostic_monthly.csv")
    f5_rows = []
    if os.path.exists(f5_csv_path):
        f5_df = pd.read_csv(f5_csv_path)
        for _, r_f5 in f5_df.iterrows():
            f5_rows.append(
                f"| {int(r_f5['month']):02d} | {float(r_f5['ol_pct_under_wlt']):.2f}% | {float(r_f5['smap_pct_under_wlt']):.2f}% |"
            )
    f5_table = "\n".join(f5_rows)

    report_md = f"""# Rapport de Synthèse Phase 2a : Évaluation des Références (Baselines B0, B1, B2)

**Date d'exécution** : {time.strftime('%Y-%m-%d %H:%M:%S')}  
**Protocole temporel** :
- **Apprentissage (train)** : 2016–2018 (`split == 'train'`).
- **Évaluation (val)** : Année 2019 (`split == 'val'`).
- **Test scellé (test)** : Année 2020 reste scellée et non consultée.

> [!WARNING]
> **Information méthodologique sur la référence B2 (fichiers CDF LDT)** :
> L'examen des fichiers de configuration officiels LDT (`ldt.config.cdf.model.NorthMor` et `ldt.config.cdf.obs.NorthMor`) indique :
> - Date de début : `2016-01-01 00:00:00`
> - Date de fin : `2021-01-01 00:00:00`
> Les courbes CDF officielles du modèle et de SMAP intègrent les années 2019 et 2020 dans leur calcul. La référence B2 incorpore ainsi de l'information issue des années d'évaluation et de test, contrairement à la référence B1 dont les paramètres sont calculés uniquement sur 2016–2018. Aucun recalcul des fichiers CDF LDT n'a été effectué sans accord préalable.

---

## 1. Description des Références

1. **B0 (NoCDF - Pas de correction)** :
   - Assimilation de l'observation brute : $\\text{{smap}}^*_0 = \\text{{smap_sm}}$, avec $d_{{\\text{{res}}, 0}} = \\text{{smap_sm}} - \\text{{ssm_ol_daily}}$.
2. **B1 (Climatologie pixel-mois)** :
   - Moyenne de $d$ par pixel et par mois, $\\hat{{c}}(p, m)$, calculée sur la période 2016–2018 pour toutes les observations passant les règles de contrôle qualité hors flétrissement.
   - Sur 2019 : $\\text{{smap}}^*_1 = \\text{{smap_sm}} - \\hat{{c}}(p, m)$ et $d_{{\\text{{res}}, 1}} = d - \\hat{{c}}(p, m)$.
3. **B2 (CDF matching v2 hors ligne)** :
   - Application de la transformation `LIS_rescale_with_CDF_matching` avec les fichiers LDT `noahmp_model_cdf_NorthMor.nc` et `smap_obs_cdf_NorthMor.nc`.
   - Sur 2019 : $\\text{{smap}}^*_2 = \\text{{CDF_rescale}}(\\text{{smap_sm}}, p, m)$ et $d_{{\\text{{res}}, 2}} = \\text{{smap}}^*_2 - \\text{{ssm_ol_daily}}$.

---

## 2. Vérification Numérique du CDF Matching contre DAOBS LIS (Record 2)

Confrontation de la fonction hors ligne avec l'enregistrement 2 des DAOBS sur 5 dates échantillonnées :

| Date et Heure DAOBS | Nombre de points comparés | Biais absolu moyen | Différence absolue max | Concordance exacte (< 1e-5) |
|:---|:---:|:---:|:---:|:---:|
{daobs_table}

---

## 3. Données Non Mises à l'Échelle par B2 (IQR et Hors Plage)

Nombre et proportion d'observations du socle hors flétrissement ({n_base_qc:,} observations) pour lesquelles la transformation CDF B2 ne renvoie pas de valeur numérique valide (critère IQR < 0.05 de la plage dynamique, ou valeur en dehors de l'intervalle [0.01, 1.0]) :

| Période / Saison | Observations socle (hors flétrissement) | B2 valeurs valides | B2 part valide (%) | B2 non transformé (NaN / IQR) | B2 part non transformée (%) |
|:---|:---:|:---:|:---:|:---:|:---:|
{b2_miss_table}

---

## 4. Population d'Évaluation et Observations Assimilables Après Correction

Pour chaque référence, le seuil de flétrissement ($\\text{{smap}}^* < \\theta_{{\\text{{wlt}}}} + 0.02$) est recalculé sur l'observation issue de la correction. Une observation est assimilable si elle satisfait les contrôles hors flétrissement, possède une valeur corrigée valide (non-NaN, comprise entre 0.01 et 1.0) et respecte $\\text{{smap}}^* \\ge \\theta_{{\\text{{wlt}}}} + 0.02$ :

{pop_table}

---

## 5. Performances Globales sur l'Année 2019

Calculées sur l'échantillon assimilable propre à chaque méthode, ainsi que sur l'échantillon commun aux trois références ({glob['n_common']:,} observations) :

{glob_metrics_table}

- **B0** : Biais de {glob['bias_b0_spec']:+.4f} m³/m³ (spécifique) et {glob['bias_b0_com']:+.4f} m³/m³ (commun) ; RMSE de {glob['rmse_b0_spec']:.4f} m³/m³ (spécifique) et {glob['rmse_b0_com']:.4f} m³/m³ (commun).
- **B1** : Biais de {glob['bias_b1_spec']:+.4f} m³/m³ (spécifique) et {glob['bias_b1_com']:+.4f} m³/m³ (commun) ; RMSE de {glob['rmse_b1_spec']:.4f} m³/m³ (spécifique) et {glob['rmse_b1_com']:.4f} m³/m³ (commun).
- **B2** : Biais de {glob['bias_b2_spec']:+.4f} m³/m³ (spécifique) et {glob['bias_b2_com']:+.4f} m³/m³ (commun) ; RMSE de {glob['rmse_b2_spec']:.4f} m³/m³ (spécifique) et {glob['rmse_b2_com']:.4f} m³/m³ (commun).

---

## 6. Corrélation d'Anomalies avec SMAP Brut

Corrélation de Pearson entre les anomalies temporelles de l'observation corrigée et celles de SMAP brut (anomalie = observation moins moyenne locale du pixel pour le mois considéré) :

{anom_table}

---

## 7. Résultats Stratifiés par Saison (Année 2019)

| Saison | RMSE B0 (spécifique) | RMSE B1 (spécifique) | RMSE B2 (spécifique) | N Obs Commun | RMSE B0 (commun) | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{season_table}

---

## 8. Résultats Stratifiés par Texture Dominante (Top 5 Classes)

| Texture | N Obs B0 | RMSE B0 | N Obs B1 | RMSE B1 | N Obs B2 | RMSE B2 | N Obs Commun | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{tex_table}

---

## 9. Résultats Stratifiés par Classe d'Irrigation GMIA

| Classe GMIA | N Obs B0 | RMSE B0 | N Obs B1 | RMSE B1 | N Obs B2 | RMSE B2 | N Obs Commun | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{gmia_table}

---

## 10. Contraste du Résidu entre GMIA >30 % et GMIA 0 % à Texture Fixée

Écart $\\Delta d_{{\\text{{res}}}} = \\bar{{d}}_{{\\text{{res}}}}(\\text{{GMIA}} > 30\\%) - \\bar{{d}}_{{\\text{{res}}}}(\\text{{GMIA}} = 0\\%)$ calculé à texture fixée (Clay Loam et Loam), avec intervalle de confiance à 95 % par rééchantillonnage bootstrap sur blocs spatiaux de 1° ($B=500$) et nombre de pixels distincts :

| Texture | Saison | Échantillon | Contraste B0 [IC 95%] (px >30%/0%) | Contraste B1 [IC 95%] (px >30%/0%) | Contraste B2 [IC 95%] (px >30%/0%) |
|:---|:---|:---|:---:|:---:|:---:|
{contrast_table}

---

## 11. Diagnostic du Dessèchement Estival (Figure F5)

Part des valeurs sous le seuil $\\theta_{{\\text{{wlt}}}} + 0.02$ calculée par mois pour le modèle Noah-MP Open Loop (`ssm_ol_daily`) et pour SMAP brut (`smap_sm`) sur la série 2016–2020 :

| Mois | Noah-MP Open Loop (% sous seuil) | SMAP brut (% sous seuil) |
|:---:|:---:|:---:|
{f5_table}

- Figure associée : [F5_wilt_diagnostic_ol_vs_smap.png](file:///home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures/F5_wilt_diagnostic_ol_vs_smap.png) (version anglaise : [F5_wilt_diagnostic_ol_vs_smap_en.png](file:///home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures/en/F5_wilt_diagnostic_ol_vs_smap_en.png)).
- Tableau de données mensuel : [figures/descriptive/F5_wilt_diagnostic_monthly.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F5_wilt_diagnostic_monthly.csv).
- Données spatiales JJA : [figures/descriptive/F5_wilt_diagnostic_jja_map.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F5_wilt_diagnostic_jja_map.csv).

---

## 12. Effectifs de la Figure F3 (Nombre d'Observations et Pixels Distincts)

Extrait des effectifs sous la figure F3 issus de [figures/descriptive/F3_d_by_gmia_texture_season.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F3_d_by_gmia_texture_season.csv) :

| Texture | Période | Classe GMIA | Nombre d'observations | Pixels distincts | Moyenne ± demi-IC 95% |
|:---|:---|:---:|:---:|:---:|:---:|
{f3_table}
"""

    report_path = os.path.join(REPORTS_DIR, "phase2a_baselines_report.md")
    with open(report_path, "w") as f:
        f.write(report_md)
    print(f"Generated synthesis report: {report_path}")


def main():
    print("=== STARTING PHASE 2a REFINED BASELINES BENCHMARK ===")
    t_start = time.time()

    # Load static masks and predictors
    with nc.Dataset(STATIC_FILE) as ds:
        landmask = ds.variables['landmask'][:]
        soil_wlt = ds.variables['soil_wlt'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]

    # 1. Load or build CDF lookup table
    lut = load_or_build_cdf_lut(landmask)

    # 2. Verify CDF matching against DAOBS Record 2
    daobs_verif = verify_cdf_matching_against_daobs(lut, landmask)

    # 3. Learn B1 climatology on 2016-2018 train split with non-wilt QC
    climatology, pixel_fallback, global_fallback = learn_b1_climatology()

    # 4. Evaluate B0, B1, B2 on 2019 validation set starting from non-wilt QC
    df_2019, n_total_2019, n_base_qc = evaluate_baselines_on_2019(
        lut, climatology, pixel_fallback, global_fallback, landmask, soil_wlt, dom_soil
    )

    # 5. Compute B2 missing/unrescaled stats
    b2_missing_df = compute_b2_missing_stats(df_2019)

    # 6. Compute anomaly correlations
    anom_corr = compute_anomaly_correlations(df_2019)

    # 7. Compute GMIA contrast (>30% vs 0%) with 1° block bootstrap
    contrast_df = compute_gmia_contrast_with_bootstrap(df_2019)

    # 8. Compute comprehensive stratified metrics (specific and common)
    metrics_df = compute_comprehensive_strata_metrics(df_2019, n_total_2019, n_base_qc)

    # 9. Generate report and CSVs
    generate_refined_synthesis_report(
        metrics_df, b2_missing_df, contrast_df, anom_corr, daobs_verif, n_total_2019, n_base_qc
    )

    print(f"=== PHASE 2a COMPLETED SUCCESSFULLY IN {time.time()-t_start:.1f}s ===")


if __name__ == "__main__":
    main()
