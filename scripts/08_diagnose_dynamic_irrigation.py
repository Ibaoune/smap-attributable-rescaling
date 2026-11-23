#!/usr/bin/env python3
"""
scripts/08_diagnose_dynamic_irrigation.py

PHASE 2a-bis: Diagnostic of Dynamic Irrigation Signal (2016-2019, 2020 strictly sealed).
Protocol:
1. Identify wetting events without precipitation:
   - Successive pairs of observations from the same pass (AM or PM) at the same pixel,
     spaced by at most 3 days (1 <= delta_t <= 3).
   - Observed increase: delta_SM > threshold_sm (default: 0.03 m3/m3).
   - Cumulative IMERG precipitation between the two observations < threshold_rain (default: 1.0 mm).
   - Sensitivity analysis on 2 other threshold combinations (lax: 0.02 m3/m3 / 1.0 mm, strict: 0.04 m3/m3 / 0.5 mm).
2. Rigorous evaluation replacing the uninformative zero-OL control:
   (a) Difference-in-Differences (DiD):
       [freq(>30%) - freq(0%)]_JJA - [freq(>30%) - freq(0%)]_DJF
       at fixed texture, with 95% CI via 1° spatial block bootstrap and distinct pixel + block counts.
   (b) Continuous GMIA trend test:
       Linear trend slope beta (% frequency per 100% GMIA) and Pearson correlation r,
       with 95% CI via 1° spatial block bootstrap across seasons and textures.
3. Outputs:
   - reports/phase2a_irrigation_diagnostic.csv
   - reports/phase2a_irrigation_did.csv
   - reports/phase2a_irrigation_trend.csv
   - reports/phase2a_irrigation_diagnostic_sensitivity.csv
   - reports/phase2a_irrigation_diagnostic_report.md (purely computed facts, zero hardcoded numbers)
   - figures/descriptive/F6_dynamic_irrigation_wetting_events.png (and en/ version)
"""

import os
import time
import yaml
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Paths (Relative to repository root)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
PRECIP_DIR = os.path.join(BASE_DIR, "data", "interim", "noahmp", "forcing", "precipitation")
CONFIG_FILE = os.path.join(BASE_DIR, "configs", "irrigation_diagnostic.yaml")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
FIG_DIR = os.path.join(BASE_DIR, "figures", "descriptive")
FIG_EN_DIR = os.path.join(FIG_DIR, "en")

TEXTURE_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def load_config():
    with open(CONFIG_FILE, 'r') as f:
        cfg = yaml.safe_load(f)
    return cfg['thresholds']


def load_continuous_imerg():
    print("Loading continuous IMERG precipitation for 2016-2019...")
    t0 = time.time()
    imerg_list = []
    for y in [2016, 2017, 2018, 2019]:
        fpath = os.path.join(PRECIP_DIR, f"precipitation_{y}.nc")
        with xr.open_dataset(fpath) as ds:
            imerg_list.append(ds['precipitation'].values)
    imerg = np.concatenate(imerg_list, axis=0)  # (1461, 121, 181)
    imerg_flat = imerg.reshape(len(imerg), 121 * 181)
    print(f"Loaded IMERG array: shape {imerg.shape} in {time.time()-t0:.2f}s")
    return imerg_flat


def load_and_pair_observations(imerg_flat):
    print("Loading 2016-2019 observations and pairing consecutive passes...")
    t0 = time.time()
    
    with nc.Dataset(STATIC_FILE) as ds:
        dom_soil = ds.variables['dominant_soiltype'][:]
    dom_soil_flat = dom_soil.flatten()
    
    cols = [
        'pixel_id', 'date', 'month', 'season', 'pass', 'smap_sm',
        'irrigfrac_gmia', 'gmia_class', 'spatial_block_id',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    
    dfs = []
    for y in [2016, 2017, 2018, 2019]:
        p = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(p, columns=cols)
        qc_nw = (
            (df_y['qc_lis_rain'] == 0) &
            (df_y['qc_lis_frozen'] == 0) &
            (df_y['qc_lis_forest'] == 0) &
            (df_y['qc_lis_snow'] == 0) &
            (df_y['qc_lis_high_veg'] == 0) &
            (df_y['qc_lis_sat_limit'] == 0) &
            (df_y['qc_lis_vwc_gt_5'] == 0) &
            (~df_y['smap_sm'].isna()) &
            (df_y['smap_sm'] >= 0.01) &
            (df_y['smap_sm'] <= 1.0)
        )
        dfs.append(df_y[qc_nw])
        
    full_df = pd.concat(dfs, ignore_index=True)
    full_df['date_dt'] = pd.to_datetime(full_df['date'])
    date0 = pd.to_datetime('2016-01-01')
    full_df['day_idx'] = (full_df['date_dt'] - date0).dt.days
    
    full_df = full_df.sort_values(['pixel_id', 'pass', 'day_idx']).reset_index(drop=True)
    
    same_group = (full_df['pixel_id'] == full_df['pixel_id'].shift(1)) & (full_df['pass'] == full_df['pass'].shift(1))
    d_days = full_df['day_idx'] - full_df['day_idx'].shift(1)
    valid_pair = same_group & (d_days >= 1) & (d_days <= 3)
    
    t2 = full_df['day_idx'].values
    t1 = full_df['day_idx'].shift(1).fillna(-999).astype(int).values
    p_ids = full_df['pixel_id'].values
    dt = d_days.fillna(-999).astype(int).values
    
    rain_cum = np.zeros(len(full_df), dtype=np.float32)
    for d in [1, 2, 3]:
        idx_d = np.where(valid_pair.values & (dt == d))[0]
        if len(idx_d) == 0:
            continue
        acc = np.zeros(len(idx_d), dtype=np.float32)
        for step in range(d):
            acc += imerg_flat[t1[idx_d] + step, p_ids[idx_d]]
        rain_cum[idx_d] = acc + imerg_flat[t2[idx_d], p_ids[idx_d]]
        
    full_df['valid_pair'] = valid_pair
    full_df['delta_days'] = dt
    full_df['rain_cum_imerg'] = rain_cum
    full_df['delta_smap'] = full_df['smap_sm'] - full_df['smap_sm'].shift(1)
    
    textures = [TEXTURE_NAMES.get(dom_soil_flat[pix], "Other") for pix in full_df['pixel_id'].values]
    full_df['texture_name'] = textures
    
    pairs_df = full_df[full_df['valid_pair']].copy().reset_index(drop=True)
    print(f"Extracted {len(pairs_df):,} valid pairs in {time.time()-t0:.2f}s")
    return pairs_df


def compute_stratum_bootstrap_rate(events, totals, B=500):
    """Vectorized 1° spatial block bootstrap for event frequency rate (%)"""
    n_blocks = len(events)
    if n_blocks == 0 or totals.sum() == 0:
        return 0.0, np.nan, np.nan
    base_rate = float(events.sum() / totals.sum() * 100.0)
    if n_blocks < 2:
        return base_rate, base_rate, base_rate
        
    idx = np.random.choice(n_blocks, size=(B, n_blocks), replace=True)
    boot_ev = events[idx].sum(axis=1)
    boot_tot = totals[idx].sum(axis=1)
    valid = boot_tot > 0
    if valid.sum() < 50:
        return base_rate, np.nan, np.nan
    rates = boot_ev[valid] / boot_tot[valid] * 100.0
    ci_low, ci_high = np.percentile(rates, [2.5, 97.5])
    return base_rate, float(ci_low), float(ci_high)


def evaluate_diagnostic(pairs_df, sm_thresh, rain_thresh):
    """Evaluate wetting without rain frequencies across strata."""
    pairs_df['is_wet_smap'] = (pairs_df['delta_smap'] > sm_thresh) & (pairs_df['rain_cum_imerg'] < rain_thresh)
    
    rows = []
    np.random.seed(42)
    
    textures = ["Clay Loam", "Loam", "Sandy Clay Loam", "Sandy Loam", "Clay", "All"]
    gmia_classes = ["0%", "0-10%", "10-30%", ">30%", "All"]
    scopes = ["All", "JJA", "DJF", "MAM", "SON"] + [str(m) for m in range(1, 13)]
    
    eval_groups = []
    # 1. Overall
    eval_groups.append(("All", "All", "All"))
    # 2. GMIA classes overall
    for gmia in ["0%", "0-10%", "10-30%", ">30%"]:
        eval_groups.append((gmia, "All", "All"))
        eval_groups.append((gmia, "All", "JJA"))
        for m in range(1, 13):
            eval_groups.append((gmia, "All", str(m)))
    # 3. Clay Loam & Loam by GMIA
    for tex in ["Clay Loam", "Loam"]:
        for gmia in ["0%", "0-10%", "10-30%", ">30%", "All"]:
            eval_groups.append((gmia, tex, "All"))
            eval_groups.append((gmia, tex, "JJA"))
            eval_groups.append((gmia, tex, "DJF"))
            for m in range(1, 13):
                eval_groups.append((gmia, tex, str(m)))
                
    seen = set()
    for gmia, tex, scope in eval_groups:
        key = (gmia, tex, scope)
        if key in seen:
            continue
        seen.add(key)
        
        mask = np.ones(len(pairs_df), dtype=bool)
        if gmia != "All":
            mask &= (pairs_df['gmia_class'] == gmia)
        if tex != "All":
            mask &= (pairs_df['texture_name'] == tex)
        if scope in ["JJA", "DJF", "MAM", "SON"]:
            mask &= (pairs_df['season'] == scope)
        elif scope != "All":
            mask &= (pairs_df['month'] == int(scope))
            
        sub = pairs_df[mask]
        n_pairs = len(sub)
        n_px = int(sub['pixel_id'].nunique())
        n_blocks = int(sub['spatial_block_id'].nunique())
        
        if n_pairs == 0:
            continue
            
        block_grp = sub.groupby('spatial_block_id')['is_wet_smap'].agg(['sum', 'count'])
        ev_smap = block_grp['sum'].values
        tot = block_grp['count'].values
        
        smap_rate, smap_low, smap_high = compute_stratum_bootstrap_rate(ev_smap, tot, B=500)
        if n_blocks < 5:
            smap_low, smap_high = np.nan, np.nan
            status = 'non interprétable (< 5 blocs)'
        else:
            status = 'interprétable'
        
        rows.append({
            'gmia_class': gmia,
            'texture': tex,
            'temporal_scope': scope,
            'n_pairs': n_pairs,
            'n_pixels': n_px,
            'n_blocks': n_blocks,
            'smap_wet_events': int(ev_smap.sum()),
            'smap_freq_pct': smap_rate,
            'smap_ci_low': smap_low,
            'smap_ci_high': smap_high,
            'status': status
        })
        
    return pd.DataFrame(rows)


def compute_did_with_bootstrap(pairs_df, sm_thresh, rain_thresh, B=500):
    """Compute Difference-in-Differences (DiD) of wetting event frequency:
       [freq(>30%) - freq(0%)]_JJA - [freq(>30%) - freq(0%)]_DJF
       at fixed texture, with 1° block bootstrap 95% CIs.
    """
    print("Computing Difference-in-Differences (DiD) with 1° block bootstrap...")
    np.random.seed(42)
    pairs_df['is_wet_smap'] = ((pairs_df['delta_smap'] > sm_thresh) & (pairs_df['rain_cum_imerg'] < rain_thresh)).astype(int)
    
    target_textures = ["Clay Loam", "Loam", "Sandy Clay Loam", "Sandy Loam", "Clay", "All"]
    rows = []
    
    for tex in target_textures:
        tex_mask = (pairs_df['texture_name'] == tex) if tex != "All" else np.ones(len(pairs_df), dtype=bool)
        sub = pairs_df[tex_mask & pairs_df['season'].isin(['JJA', 'DJF']) & pairs_df['gmia_class'].isin(['0%', '>30%'])]
        
        if len(sub) == 0:
            continue
            
        # Distinct blocks and pixels
        blocks = np.sort(sub['spatial_block_id'].unique())
        n_b = len(blocks)
        if n_b < 2:
            continue
            
        b_map = {b: i for i, b in enumerate(blocks)}
        b_indices = np.array([b_map[b] for b in sub['spatial_block_id'].values])
        
        is_jja = (sub['season'].values == 'JJA')
        is_djf = (sub['season'].values == 'DJF')
        is_gt30 = (sub['gmia_class'].values == '>30%')
        is_0 = (sub['gmia_class'].values == '0%')
        is_wet = sub['is_wet_smap'].values
        
        # Accumulate sums per block for each combination
        # (JJA >30%, JJA 0%, DJF >30%, DJF 0%)
        ev_jja_gt30 = np.bincount(b_indices[is_jja & is_gt30 & (is_wet == 1)], minlength=n_b)
        tot_jja_gt30 = np.bincount(b_indices[is_jja & is_gt30], minlength=n_b)
        
        ev_jja_0 = np.bincount(b_indices[is_jja & is_0 & (is_wet == 1)], minlength=n_b)
        tot_jja_0 = np.bincount(b_indices[is_jja & is_0], minlength=n_b)
        
        ev_djf_gt30 = np.bincount(b_indices[is_djf & is_gt30 & (is_wet == 1)], minlength=n_b)
        tot_djf_gt30 = np.bincount(b_indices[is_djf & is_gt30], minlength=n_b)
        
        ev_djf_0 = np.bincount(b_indices[is_djf & is_0 & (is_wet == 1)], minlength=n_b)
        tot_djf_0 = np.bincount(b_indices[is_djf & is_0], minlength=n_b)
        
        # Base rates
        rate_jja_gt30 = ev_jja_gt30.sum() / tot_jja_gt30.sum() * 100.0 if tot_jja_gt30.sum() > 0 else np.nan
        rate_jja_0 = ev_jja_0.sum() / tot_jja_0.sum() * 100.0 if tot_jja_0.sum() > 0 else np.nan
        diff_jja = rate_jja_gt30 - rate_jja_0
        
        rate_djf_gt30 = ev_djf_gt30.sum() / tot_djf_gt30.sum() * 100.0 if tot_djf_gt30.sum() > 0 else np.nan
        rate_djf_0 = ev_djf_0.sum() / tot_djf_0.sum() * 100.0 if tot_djf_0.sum() > 0 else np.nan
        diff_djf = rate_djf_gt30 - rate_djf_0
        
        did = diff_jja - diff_djf
        
        # Block Bootstrap
        idx_matrix = np.random.choice(n_b, size=(B, n_b), replace=True)
        
        tot_jja_gt30_b = tot_jja_gt30[idx_matrix].sum(axis=1)
        tot_jja_0_b = tot_jja_0[idx_matrix].sum(axis=1)
        tot_djf_gt30_b = tot_djf_gt30[idx_matrix].sum(axis=1)
        tot_djf_0_b = tot_djf_0[idx_matrix].sum(axis=1)
        
        valid = (tot_jja_gt30_b > 0) & (tot_jja_0_b > 0) & (tot_djf_gt30_b > 0) & (tot_djf_0_b > 0)
        if valid.sum() > 50:
            diff_jja_b = (ev_jja_gt30[idx_matrix].sum(axis=1)[valid] / tot_jja_gt30_b[valid] -
                          ev_jja_0[idx_matrix].sum(axis=1)[valid] / tot_jja_0_b[valid]) * 100.0
            diff_djf_b = (ev_djf_gt30[idx_matrix].sum(axis=1)[valid] / tot_djf_gt30_b[valid] -
                          ev_djf_0[idx_matrix].sum(axis=1)[valid] / tot_djf_0_b[valid]) * 100.0
            did_b = diff_jja_b - diff_djf_b
            
            jja_low, jja_high = float(np.percentile(diff_jja_b, 2.5)), float(np.percentile(diff_jja_b, 97.5))
            djf_low, djf_high = float(np.percentile(diff_djf_b, 2.5)), float(np.percentile(diff_djf_b, 97.5))
            did_low, did_high = float(np.percentile(did_b, 2.5)), float(np.percentile(did_b, 97.5))
        else:
            jja_low, jja_high = np.nan, np.nan
            djf_low, djf_high = np.nan, np.nan
            did_low, did_high = np.nan, np.nan
            
        sub_gt30 = sub[sub['gmia_class'] == '>30%']
        sub_0 = sub[sub['gmia_class'] == '0%']
        n_b_gt30 = int(sub_gt30['spatial_block_id'].nunique()) if len(sub_gt30) > 0 else 0
        n_b_0 = int(sub_0['spatial_block_id'].nunique()) if len(sub_0) > 0 else 0
        
        if n_b_gt30 < 5 or n_b_0 < 5 or n_b < 5:
            jja_low, jja_high = np.nan, np.nan
            djf_low, djf_high = np.nan, np.nan
            did_low, did_high = np.nan, np.nan
            status = 'non interprétable (< 5 blocs)'
        else:
            status = 'interprétable'
        
        rows.append({
            'texture': tex,
            'n_pairs_jja_gt30': int(tot_jja_gt30.sum()),
            'n_pairs_jja_0': int(tot_jja_0.sum()),
            'n_pairs_djf_gt30': int(tot_djf_gt30.sum()),
            'n_pairs_djf_0': int(tot_djf_0.sum()),
            'n_px_gt30': int(sub_gt30['pixel_id'].nunique()) if len(sub_gt30) > 0 else 0,
            'n_px_0': int(sub_0['pixel_id'].nunique()) if len(sub_0) > 0 else 0,
            'n_blocks_gt30': n_b_gt30,
            'n_blocks_0': n_b_0,
            'n_blocks_total': n_b,
            'rate_jja_gt30': rate_jja_gt30,
            'rate_jja_0': rate_jja_0,
            'diff_jja': diff_jja,
            'diff_jja_ci_low': jja_low,
            'diff_jja_ci_high': jja_high,
            'rate_djf_gt30': rate_djf_gt30,
            'rate_djf_0': rate_djf_0,
            'diff_djf': diff_djf,
            'diff_djf_ci_low': djf_low,
            'diff_djf_ci_high': djf_high,
            'did': did,
            'did_ci_low': did_low,
            'did_ci_high': did_high,
            'status': status
        })
        
    return pd.DataFrame(rows)


def compute_continuous_trend_with_bootstrap(pairs_df, sm_thresh, rain_thresh, B=500):
    """Compute continuous GMIA trend test of wetting event probability (expressed in %)
       against continuous GMIA fraction [0, 1] across strata with 1° block bootstrap 95% CIs.
    """
    print("Computing continuous GMIA trend tests with 1° block bootstrap...")
    np.random.seed(42)
    pairs_df['is_wet_smap'] = ((pairs_df['delta_smap'] > sm_thresh) & (pairs_df['rain_cum_imerg'] < rain_thresh)).astype(int)
    
    strata = [
        ("Clay Loam", "JJA"),
        ("Clay Loam", "DJF"),
        ("Clay Loam", "All"),
        ("Loam", "JJA"),
        ("Loam", "DJF"),
        ("Loam", "All"),
        ("All", "JJA"),
        ("All", "DJF"),
        ("All", "All")
    ]
    
    rows = []
    for tex, s in strata:
        mask = np.ones(len(pairs_df), dtype=bool)
        if tex != "All":
            mask &= (pairs_df['texture_name'] == tex)
        if s != "All":
            mask &= (pairs_df['season'] == s)
            
        sub = pairs_df[mask]
        n_pairs = len(sub)
        if n_pairs == 0:
            continue
            
        blocks = np.sort(sub['spatial_block_id'].unique())
        n_b = len(blocks)
        if n_b < 2:
            continue
            
        b_map = {b: i for i, b in enumerate(blocks)}
        b_indices = np.array([b_map[b] for b in sub['spatial_block_id'].values])
        
        x = sub['irrigfrac_gmia'].values.astype(np.float64)
        y = sub['is_wet_smap'].values.astype(np.float64) * 100.0  # Percentage
        
        # Block-level sums
        s0_b = np.bincount(b_indices, minlength=n_b).astype(np.float64)
        sx_b = np.bincount(b_indices, weights=x, minlength=n_b)
        sy_b = np.bincount(b_indices, weights=y, minlength=n_b)
        sxx_b = np.bincount(b_indices, weights=x*x, minlength=n_b)
        syy_b = np.bincount(b_indices, weights=y*y, minlength=n_b)
        sxy_b = np.bincount(b_indices, weights=x*y, minlength=n_b)
        
        s0 = s0_b.sum()
        sx = sx_b.sum()
        sy = sy_b.sum()
        sxx = sxx_b.sum()
        syy = syy_b.sum()
        sxy = sxy_b.sum()
        
        denom_x = s0 * sxx - sx**2
        denom_y = s0 * syy - sy**2
        
        if denom_x > 1e-12:
            beta = (s0 * sxy - sx * sy) / denom_x
            r = (s0 * sxy - sx * sy) / np.sqrt(denom_x * denom_y) if denom_y > 1e-12 else 0.0
        else:
            beta = 0.0
            r = 0.0
            
        # Bootstrap
        idx_matrix = np.random.choice(n_b, size=(B, n_b), replace=True)
        s0_m = s0_b[idx_matrix].sum(axis=1)
        sx_m = sx_b[idx_matrix].sum(axis=1)
        sy_m = sy_b[idx_matrix].sum(axis=1)
        sxx_m = sxx_b[idx_matrix].sum(axis=1)
        syy_m = syy_b[idx_matrix].sum(axis=1)
        sxy_m = sxy_b[idx_matrix].sum(axis=1)
        
        den_x_m = s0_m * sxx_m - sx_m**2
        den_y_m = s0_m * syy_m - sy_m**2
        valid = (den_x_m > 1e-12) & (s0_m > 0)
        
        if valid.sum() > 50:
            beta_m = (s0_m[valid] * sxy_m[valid] - sx_m[valid] * sy_m[valid]) / den_x_m[valid]
            r_m = (s0_m[valid] * sxy_m[valid] - sx_m[valid] * sy_m[valid]) / np.sqrt(den_x_m[valid] * den_y_m[valid])
            b_low, b_high = float(np.percentile(beta_m, 2.5)), float(np.percentile(beta_m, 97.5))
            r_low, r_high = float(np.percentile(r_m, 2.5)), float(np.percentile(r_m, 97.5))
        else:
            b_low, b_high = np.nan, np.nan
            r_low, r_high = np.nan, np.nan
            
        if n_b < 5:
            b_low, b_high = np.nan, np.nan
            r_low, r_high = np.nan, np.nan
            status = 'non interprétable (< 5 blocs)'
        else:
            status = 'interprétable'

        rows.append({
            'stratum_texture': tex,
            'stratum_season': s,
            'n_pairs': n_pairs,
            'n_pixels': int(sub['pixel_id'].nunique()),
            'n_blocks': n_b,
            'mean_gmia_pct': float(x.mean() * 100.0),
            'mean_wet_freq_pct': float(y.mean()),
            'beta_slope': float(beta),
            'beta_ci_low': b_low,
            'beta_ci_high': b_high,
            'pearson_r': float(r),
            'r_ci_low': r_low,
            'r_ci_high': r_high,
            'status': status
        })
        
    return pd.DataFrame(rows)


def run_sensitivity(pairs_df, thresh_cfg):
    """Run sensitivity analysis across configured threshold combinations."""
    print("Running sensitivity analysis across threshold configurations...")
    records = []
    
    for cfg_key in ['default', 'sensitivity_lax', 'sensitivity_strict']:
        cfg = thresh_cfg[cfg_key]
        sm_th = cfg['delta_sm_min']
        rain_th = cfg['rain_imerg_max']
        
        pairs_df['is_wet'] = (pairs_df['delta_smap'] > sm_th) & (pairs_df['rain_cum_imerg'] < rain_th)
        
        # Overall
        ov_smap = pairs_df['is_wet'].mean() * 100.0
        
        # Clay Loam JJA
        cl_jja = pairs_df[(pairs_df['texture_name'] == 'Clay Loam') & (pairs_df['season'] == 'JJA')]
        cl_gt30 = cl_jja[cl_jja['gmia_class'] == '>30%']['is_wet'].mean() * 100.0
        cl_0 = cl_jja[cl_jja['gmia_class'] == '0%']['is_wet'].mean() * 100.0
        diff_cl = cl_gt30 - cl_0
        
        # Loam JJA
        loam_jja = pairs_df[(pairs_df['texture_name'] == 'Loam') & (pairs_df['season'] == 'JJA')]
        loam_gt30 = loam_jja[loam_jja['gmia_class'] == '>30%']['is_wet'].mean() * 100.0
        loam_0 = loam_jja[loam_jja['gmia_class'] == '0%']['is_wet'].mean() * 100.0
        diff_loam = loam_gt30 - loam_0
        
        records.append({
            'threshold_config': cfg_key,
            'delta_sm_min': sm_th,
            'rain_imerg_max': rain_th,
            'overall_smap_freq_pct': ov_smap,
            'clay_loam_jja_gt30_pct': cl_gt30,
            'clay_loam_jja_0_pct': cl_0,
            'clay_loam_jja_contrast_pct': diff_cl,
            'loam_jja_gt30_pct': loam_gt30,
            'loam_jja_0_pct': loam_0,
            'loam_jja_contrast_pct': diff_loam
        })
        
    return pd.DataFrame(records)


def plot_f6_figure(df_default, df_did, df_trend, df_sens, pairs_df, thresh_cfg, out_dir, lang='fr'):
    """Generate publication figure F6 showing monthly cycle, DiD comparison, continuous trend, and sensitivity."""
    os.makedirs(out_dir, exist_ok=True)
    fname = "F6_dynamic_irrigation_wetting_events.png" if lang == 'fr' else "F6_dynamic_irrigation_wetting_events_en.png"
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # Panel A: Monthly cycle in Clay Loam (0% vs >30%)
    ax0 = axes[0, 0]
    months = np.arange(1, 13)
    cl_months = df_default[(df_default['texture'] == 'Clay Loam') & (df_default['temporal_scope'].isin([str(m) for m in months]))]
    
    sub_gt30 = cl_months[cl_months['gmia_class'] == '>30%'].sort_values('temporal_scope', key=lambda x: x.astype(int))
    sub_0 = cl_months[cl_months['gmia_class'] == '0%'].sort_values('temporal_scope', key=lambda x: x.astype(int))
    
    ax0.plot(months, sub_gt30['smap_freq_pct'].values, 'o-', color='#d7191c', linewidth=2, label='GMIA > 30% (Clay Loam)')
    ax0.fill_between(months, sub_gt30['smap_ci_low'].values, sub_gt30['smap_ci_high'].values, color='#d7191c', alpha=0.15)
    
    ax0.plot(months, sub_0['smap_freq_pct'].values, 's--', color='#2b83ba', linewidth=1.5, label='GMIA = 0% (Clay Loam)')
    ax0.fill_between(months, sub_0['smap_ci_low'].values, sub_0['smap_ci_high'].values, color='#2b83ba', alpha=0.15)
    
    ax0.axvspan(6, 8, color='yellow', alpha=0.15, label="Été (JJA)" if lang == 'fr' else "Summer (JJA)")
    ax0.set_xlabel("Mois" if lang == 'fr' else "Month", fontsize=11, fontweight='bold')
    ax0.set_ylabel("Fréquence d'événements sans pluie (%)" if lang == 'fr' else "Wetting event frequency (%)", fontsize=11, fontweight='bold')
    ax0.set_title("(a) Cycle mensuel en Clay Loam (>30% vs 0%)" if lang == 'fr' else "(a) Monthly cycle in Clay Loam (>30% vs 0%)", fontsize=11)
    ax0.set_xticks(months)
    ax0.grid(True, linestyle='--', alpha=0.5)
    ax0.legend(frameon=True, fontsize=9.5)
    
    # Panel B: Difference of Differences across textures
    ax1 = axes[0, 1]
    textures_did = ["Clay Loam", "Loam", "Sandy Clay Loam", "Sandy Loam", "Clay", "All"]
    sub_did = df_did[df_did['status'] == 'interprétable'].set_index('texture')
    x_idx = np.arange(len(sub_did))
    width = 0.28
    
    yerr_jja = [(sub_did['diff_jja'] - sub_did['diff_jja_ci_low']).values, (sub_did['diff_jja_ci_high'] - sub_did['diff_jja']).values]
    yerr_djf = [(sub_did['diff_djf'] - sub_did['diff_djf_ci_low']).values, (sub_did['diff_djf_ci_high'] - sub_did['diff_djf']).values]
    yerr_did = [(sub_did['did'] - sub_did['did_ci_low']).values, (sub_did['did_ci_high'] - sub_did['did']).values]
    
    ax1.bar(x_idx - width, sub_did['diff_jja'].values, width, yerr=yerr_jja, capsize=3, color='#d7191c', label=r"$\Delta_{\mathrm{JJA}}$ (>30% - 0%)")
    ax1.bar(x_idx, sub_did['diff_djf'].values, width, yerr=yerr_djf, capsize=3, color='#2b83ba', label=r"$\Delta_{\mathrm{DJF}}$ (>30% - 0%)")
    ax1.bar(x_idx + width, sub_did['did'].values, width, yerr=yerr_did, capsize=3, color='#7b3294', label=r"DiD ($\Delta_{\mathrm{JJA}} - \Delta_{\mathrm{DJF}}$)")
    
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.8)
    ax1.set_xlabel("Texture du sol" if lang == 'fr' else "Soil texture", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Différence de fréquence (%)" if lang == 'fr' else "Frequency difference (%)", fontsize=11, fontweight='bold')
    ax1.set_title("(b) Différence de différences (Été vs Hiver)" if lang == 'fr' else "(b) Difference-in-Differences (Summer vs Winter)", fontsize=11)
    ax1.set_xticks(x_idx)
    ax1.set_xticklabels(sub_did.index, rotation=20, ha='right', fontsize=9.5)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(frameon=True, fontsize=9)
    
    # Panel C: Continuous GMIA trend (JJA vs DJF in Clay Loam)
    ax2 = axes[1, 0]
    cl_pairs = pairs_df[pairs_df['texture_name'] == 'Clay Loam']
    cl_jja = cl_pairs[cl_pairs['season'] == 'JJA']
    cl_djf = cl_pairs[cl_pairs['season'] == 'DJF']
    
    # 10 GMIA bins for visualization
    bins = np.linspace(0, 1.0, 11)
    bin_centers = 0.5 * (bins[:-1] + bins[1:]) * 100.0
    
    cl_jja_binned = cl_jja.groupby(pd.cut(cl_jja['irrigfrac_gmia'], bins, include_lowest=True), observed=False)['is_wet_smap'].mean() * 100.0
    cl_djf_binned = cl_djf.groupby(pd.cut(cl_djf['irrigfrac_gmia'], bins, include_lowest=True), observed=False)['is_wet_smap'].mean() * 100.0
    
    ax2.plot(bin_centers, cl_jja_binned.values, 'o-', color='#d7191c', linewidth=2, label="Clay Loam (JJA)")
    ax2.plot(bin_centers, cl_djf_binned.values, 's--', color='#2b83ba', linewidth=1.5, label="Clay Loam (DJF)")
    
    # Trend line annotations from df_trend
    t_jja = df_trend[(df_trend['stratum_texture'] == 'Clay Loam') & (df_trend['stratum_season'] == 'JJA')].iloc[0]
    t_djf = df_trend[(df_trend['stratum_texture'] == 'Clay Loam') & (df_trend['stratum_season'] == 'DJF')].iloc[0]
    
    lbl_jja = f"Pente JJA : β = {t_jja['beta_slope']:+.2f}% / 100% GMIA" if lang == 'fr' else f"JJA slope: β = {t_jja['beta_slope']:+.2f}% / 100% GMIA"
    lbl_djf = f"Pente DJF : β = {t_djf['beta_slope']:+.2f}% / 100% GMIA" if lang == 'fr' else f"DJF slope: β = {t_djf['beta_slope']:+.2f}% / 100% GMIA"
    
    ax2.text(0.05, 0.90, lbl_jja, transform=ax2.transAxes, color='#d7191c', fontsize=9.5, fontweight='bold')
    ax2.text(0.05, 0.83, lbl_djf, transform=ax2.transAxes, color='#2b83ba', fontsize=9.5, fontweight='bold')
    
    ax2.set_xlabel("Fraction d'irrigation GMIA (%)" if lang == 'fr' else "GMIA irrigation fraction (%)", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Fréquence d'humidification (%)" if lang == 'fr' else "Wetting frequency (%)", fontsize=11, fontweight='bold')
    ax2.set_title("(c) Tendance continue en fonction de GMIA" if lang == 'fr' else "(c) Continuous trend vs GMIA fraction", fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='center right', frameon=True, fontsize=9.5)
    
    # Panel D: Sensitivity Analysis across threshold configurations
    ax3 = axes[1, 1]
    cfg_names = ["Laxe", "Défaut", "Stricte"] if lang == 'fr' else ["Lax", "Default", "Strict"]
    x_cfg = np.arange(len(df_sens))
    w_cfg = 0.35
    
    ax3.bar(x_cfg - w_cfg/2, df_sens['clay_loam_jja_contrast_pct'].values, w_cfg, color='#2b83ba', label="Clay Loam (JJA >30% - 0%)")
    ax3.bar(x_cfg + w_cfg/2, df_sens['loam_jja_contrast_pct'].values, w_cfg, color='#fdae61', label="Loam (JJA >30% - 0%)")
    
    ax3.set_xlabel("Configuration de seuils" if lang == 'fr' else "Threshold configuration", fontsize=11, fontweight='bold')
    ax3.set_ylabel("Contraste JJA >30% vs 0% (%)" if lang == 'fr' else "JJA contrast >30% vs 0% (%)", fontsize=11, fontweight='bold')
    ax3.set_title("(d) Analyse de sensibilité des seuils" if lang == 'fr' else "(d) Threshold sensitivity analysis", fontsize=11)
    ax3.set_xticks(x_cfg)
    ax3.set_xticklabels(cfg_names, fontsize=10)
    ax3.grid(True, linestyle='--', alpha=0.5)
    ax3.legend(frameon=True, fontsize=9.5)
    
    plt.tight_layout()
    fpath = os.path.join(out_dir, fname)
    plt.savefig(fpath, dpi=300)
    plt.close()
    print(f"Saved figure: {fpath}")


def generate_markdown_report(df_default, df_did, df_trend, df_sens, pairs_df, thresh_cfg):
    """Generate Phase 2a-bis Markdown report with zero hardcoded numbers and zero interpretive words."""
    print("Generating Phase 2a-bis Markdown report...")
    out_md = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic_report.md")
    
    n_pairs_tot = len(pairs_df)
    n_pixels_tot = int(pairs_df['pixel_id'].nunique())
    n_blocks_tot = int(pairs_df['spatial_block_id'].nunique())
    
    ov_default = df_default[(df_default['gmia_class'] == 'All') & 
                            (df_default['texture'] == 'All') & 
                            (df_default['temporal_scope'] == 'All')].iloc[0]
                            
    # DiD table rows
    did_rows = []
    for _, r_d in df_did.iterrows():
        if r_d['status'] != 'interprétable' or pd.isna(r_d['did_ci_low']):
            did_rows.append(
                f"| {r_d['texture']} | non interprétable | non interprétable | **non interprétable** | "
                f"{r_d['n_px_gt30']}/{r_d['n_px_0']} px | {r_d['n_blocks_gt30']}/{r_d['n_blocks_0']} blk (< 5) |"
            )
        else:
            did_rows.append(
                f"| {r_d['texture']} | {r_d['diff_jja']:+.3f}% [{r_d['diff_jja_ci_low']:+.3f}%, {r_d['diff_jja_ci_high']:+.3f}%] | "
                f"{r_d['diff_djf']:+.3f}% [{r_d['diff_djf_ci_low']:+.3f}%, {r_d['diff_djf_ci_high']:+.3f}%] | "
                f"**{r_d['did']:+.3f}%** [{r_d['did_ci_low']:+.3f}%, {r_d['did_ci_high']:+.3f}%] | "
                f"{r_d['n_px_gt30']}/{r_d['n_px_0']} px | {r_d['n_blocks_gt30']}/{r_d['n_blocks_0']} blk |"
            )
    table_did = chr(10).join(did_rows)
    
    # Continuous trend table rows
    trend_rows = []
    for _, r_t in df_trend.iterrows():
        if r_t['status'] != 'interprétable' or pd.isna(r_t['beta_ci_low']):
            trend_rows.append(
                f"| {r_t['stratum_texture']} | {r_t['stratum_season']} | {r_t['n_pairs']:,} | "
                f"non interprétable | non interprétable | {r_t['n_pixels']:,} px, {r_t['n_blocks']} blk (< 5) |"
            )
        else:
            trend_rows.append(
                f"| {r_t['stratum_texture']} | {r_t['stratum_season']} | {r_t['n_pairs']:,} | "
                f"{r_t['beta_slope']:+.3f}% [{r_t['beta_ci_low']:+.3f}%, {r_t['beta_ci_high']:+.3f}%] | "
                f"{r_t['pearson_r']:+.4f} [{r_t['r_ci_low']:+.4f}, {r_t['r_ci_high']:+.4f}] | "
                f"{r_t['n_pixels']:,} px, {r_t['n_blocks']} blk |"
            )
    table_trend = chr(10).join(trend_rows)
    
    # Monthly table for Clay Loam (>30% vs 0%)
    month_rows = []
    for m in range(1, 13):
        sub_gt30 = df_default[(df_default['texture'] == 'Clay Loam') & 
                              (df_default['temporal_scope'] == str(m)) & 
                              (df_default['gmia_class'] == '>30%')]
        sub_0 = df_default[(df_default['texture'] == 'Clay Loam') & 
                           (df_default['temporal_scope'] == str(m)) & 
                           (df_default['gmia_class'] == '0%')]
        if len(sub_gt30) > 0 and len(sub_0) > 0:
            r_gt30 = sub_gt30.iloc[0]
            r_0 = sub_0.iloc[0]
            if r_gt30['status'] != 'interprétable' or r_0['status'] != 'interprétable' or pd.isna(r_gt30['smap_ci_low']) or pd.isna(r_0['smap_ci_low']):
                month_rows.append(
                    f"| {m:02d} | non interprétable | non interprétable | "
                    f"non interprétable | {r_gt30['n_pixels']}/{r_0['n_pixels']} | {r_gt30['n_blocks']}/{r_0['n_blocks']} (< 5) |"
                )
            else:
                diff = r_gt30['smap_freq_pct'] - r_0['smap_freq_pct']
                month_rows.append(
                    f"| {m:02d} | {r_gt30['smap_freq_pct']:.3f}% [{r_gt30['smap_ci_low']:.3f}%, {r_gt30['smap_ci_high']:.3f}%] | "
                    f"{r_0['smap_freq_pct']:.3f}% [{r_0['smap_ci_low']:.3f}%, {r_0['smap_ci_high']:.3f}%] | "
                    f"{diff:+.3f}% | {r_gt30['n_pixels']}/{r_0['n_pixels']} | {r_gt30['n_blocks']}/{r_0['n_blocks']} |"
                )
    table_months = chr(10).join(month_rows)
    
    # Sensitivity summary table
    sens_rows = []
    for _, r_s in df_sens.iterrows():
        sens_rows.append(
            f"| {r_s['threshold_config']} | ΔSM > {r_s['delta_sm_min']:.2f} m³/m³, Rain < {r_s['rain_imerg_max']:.1f} mm | "
            f"{r_s['clay_loam_jja_gt30_pct']:.3f}% | {r_s['clay_loam_jja_0_pct']:.3f}% | {r_s['clay_loam_jja_contrast_pct']:+.3f}% | "
            f"{r_s['loam_jja_gt30_pct']:.3f}% | {r_s['loam_jja_0_pct']:.3f}% | {r_s['loam_jja_contrast_pct']:+.3f}% |"
        )
    table_sens = chr(10).join(sens_rows)
    
    latex_tau_sm = r"$\tau_{\mathrm{sm}}$"
    latex_tau_rain = r"$\tau_{\mathrm{rain}}$"
    latex_event = r"$\Delta \mathrm{SM} > \tau_{\mathrm{sm}}$"
    latex_rain = r"$\mathrm{Rain}_{t_1 \to t_2} < \tau_{\mathrm{rain}}$"
    latex_did_formula = r"$\Delta_{\mathrm{DiD}} = [f(>30\%) - f(0\%)]_{\mathrm{JJA}} - [f(>30\%) - f(0\%)]_{\mathrm{DJF}}$"
    latex_beta = r"$\beta$"
    latex_r = r"$r$"
    report_md = f"""# Rapport de Diagnostic Phase 2a-bis : Signal d'Irrigation Dynamique

**Période d'évaluation** : 2016-01-01 au 2019-12-31 (4 années complètes).  
**Statut de l'année 2020** : 2020 reste strictement scellée.  
**Population d'observation** : Paires d'observations SMAP successives du même passage orbital (AM avec AM, PM avec PM), espacées de 1 à 3 jours, issues des observations conformes aux règles QC hors flétrissement.  
**Effectif total analysé** : {n_pairs_tot:,} paires successives ({n_pixels_tot:,} mailles distinctes, {n_blocks_tot:,} blocs spatiaux distincts de 1°).

---

## 1. Définition Opérationnelle et Seuils

Un événement d'humidification sans précipitation est défini par :
- Un intervalle temporel de 1 à 3 jours entre deux observations du même passage au même pixel.
- Une élévation d'humidité du sol : {latex_event}.
- Un cumul de précipitations IMERG sur l'intervalle : {latex_rain}.
- Configuration par défaut : {latex_tau_sm} = {thresh_cfg["default"]["delta_sm_min"]:.2f} m³/m³, {latex_tau_rain} = {thresh_cfg["default"]["rain_imerg_max"]:.1f} mm.

---

## 2. Statistiques Globales (Configuration par Défaut)

- **Nombre total de paires de passages successifs** : {n_pairs_tot:,}
- **Événements détectés sur SMAP (ensemble du domaine et de l'année)** : {ov_default['smap_wet_events']:,} ({ov_default['smap_freq_pct']:.3f}%, IC 95% : [{ov_default['smap_ci_low']:.3f}%, {ov_default['smap_ci_high']:.3f}%])
- **Règle de fiabilité spatiale** : Conformément au protocole, toute strate présentant moins de 5 blocs spatiaux distincts de 1° est classée comme « non interprétable » (aucun IC n'est rapporté en raison de l'absence de degrés de liberté spatiaux suffisants).

---

## 3. Différence de Différences (DiD) : [Fréquence >30 % − 0 %] en Été (JJA) vs Hiver (DJF)

Contraste {latex_did_formula} avec intervalles de confiance à 95% par bootstrap par blocs de 1° (500 itérations) :

| Texture | Contraste JJA (>30% − 0%) | Contraste DJF (>30% − 0%) | Différence de Différences (DiD) | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---|:---:|:---:|:---:|:---:|:---:|
{table_did}

---

## 4. Test de Tendance Continue : Fréquence d'Humidification vs Fraction GMIA Continue

Pente de régression linéaire {latex_beta} (% de fréquence par 100% de GMIA) et corrélation de Pearson {latex_r} avec intervalles de confiance à 95% par bootstrap par blocs de 1° :

| Texture | Saison | Paires analysées | Pente linéaire {latex_beta} [IC 95%] | Corrélation {latex_r} [IC 95%] | Échantillon spatial |
|:---|:---:|:---:|:---:|:---:|:---:|
{table_trend}

---

## 5. Cycle Mensuel en Clay Loam : Contraste GMIA >30 % vs 0 %

| Mois | SMAP >30% [IC 95%] | SMAP 0% [IC 95%] | Contraste SMAP (>30% − 0%) | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---:|:---:|:---:|:---:|:---:|:---:|
{table_months}

---

## 6. Analyse de Sensibilité aux Seuils

| Configuration | Seuils appliqués | Clay Loam >30% | Clay Loam 0% | Contraste Clay Loam | Loam >30% | Loam 0% | Contraste Loam |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|
{table_sens}

---

## 7. Fichiers et Données Produites

- **Figure F6** : `figures/descriptive/F6_dynamic_irrigation_wetting_events.png` (version EN : `figures/descriptive/en/F6_dynamic_irrigation_wetting_events_en.png`)
- **Table complète par strate** : `reports/phase2a_irrigation_diagnostic.csv`
- **Table Différence de Différences (DiD)** : `reports/phase2a_irrigation_did.csv`
- **Table test de tendance continue GMIA** : `reports/phase2a_irrigation_trend.csv`
- **Table d'analyse de sensibilité** : `reports/phase2a_irrigation_diagnostic_sensitivity.csv`
"""
    with open(out_md, 'w') as f:
        f.write(report_md)
    print(f"Generated report: {out_md}")


def main():
    print("=== STARTING PHASE 2a-bis DYNAMIC IRRIGATION DIAGNOSTIC ===")
    t_start = time.time()
    
    thresh_cfg = load_config()
    def_cfg = thresh_cfg['default']
    
    imerg_flat = load_continuous_imerg()
    pairs_df = load_and_pair_observations(imerg_flat)
    
    # 1. Stratified rates
    print(f"Evaluating default thresholds: delta_sm > {def_cfg['delta_sm_min']}, rain < {def_cfg['rain_imerg_max']} mm...")
    df_default = evaluate_diagnostic(pairs_df, def_cfg['delta_sm_min'], def_cfg['rain_imerg_max'])
    csv_default = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic.csv")
    df_default.to_csv(csv_default, index=False)
    print(f"Saved: {csv_default}")
    
    # 2. Difference of Differences (DiD)
    df_did = compute_did_with_bootstrap(pairs_df, def_cfg['delta_sm_min'], def_cfg['rain_imerg_max'])
    csv_did = os.path.join(REPORTS_DIR, "phase2a_irrigation_did.csv")
    df_did.to_csv(csv_did, index=False)
    print(f"Saved: {csv_did}")
    
    # 3. Continuous Trend Test
    df_trend = compute_continuous_trend_with_bootstrap(pairs_df, def_cfg['delta_sm_min'], def_cfg['rain_imerg_max'])
    csv_trend = os.path.join(REPORTS_DIR, "phase2a_irrigation_trend.csv")
    df_trend.to_csv(csv_trend, index=False)
    print(f"Saved: {csv_trend}")
    
    # 4. Sensitivity analysis
    df_sens = run_sensitivity(pairs_df, thresh_cfg)
    csv_sens = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic_sensitivity.csv")
    df_sens.to_csv(csv_sens, index=False)
    print(f"Saved: {csv_sens}")
    
    # 5. Plot Figures F6 (FR & EN)
    print("Generating Figure F6 (FR and EN)...")
    plot_f6_figure(df_default, df_did, df_trend, df_sens, pairs_df, thresh_cfg, FIG_DIR, lang='fr')
    plot_f6_figure(df_default, df_did, df_trend, df_sens, pairs_df, thresh_cfg, FIG_EN_DIR, lang='en')
    
    # Copy to IDE artifacts
    ide_art_dir = "/home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures"
    os.makedirs(os.path.join(ide_art_dir, "en"), exist_ok=True)
    import shutil
    shutil.copy(os.path.join(FIG_DIR, "F6_dynamic_irrigation_wetting_events.png"), os.path.join(ide_art_dir, "F6_dynamic_irrigation_wetting_events.png"))
    shutil.copy(os.path.join(FIG_EN_DIR, "F6_dynamic_irrigation_wetting_events_en.png"), os.path.join(ide_art_dir, "en", "F6_dynamic_irrigation_wetting_events_en.png"))
    
    # 6. Generate Markdown Report
    generate_markdown_report(df_default, df_did, df_trend, df_sens, pairs_df, thresh_cfg)
    
    print(f"=== COMPLETED IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
