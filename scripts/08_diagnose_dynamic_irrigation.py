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
2. Frequency of these events stratified by:
   - GMIA class (0%, 0-10%, 10-30%, >30%) x Texture x Month (and JJA / All seasons).
   - 95% Confidence Intervals via 1° spatial block bootstrap.
   - Distinct pixel count AND distinct 1° block count reported for every stratum.
   - Identical calculation on Open Loop (delta_OL > threshold_sm without rain), serving as negative control.
3. Outputs:
   - reports/phase2a_irrigation_diagnostic.csv
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
    imerg = np.concatenate(imerg_list, axis=0) # (1461, 121, 181)
    imerg_flat = imerg.reshape(len(imerg), 121 * 181)
    print(f"Loaded IMERG array: shape {imerg.shape} in {time.time()-t0:.2f}s")
    return imerg_flat


def load_and_pair_observations(imerg_flat):
    print("Loading 2016-2019 observations and pairing consecutive passes...")
    t0 = time.time()
    
    # Load texture from static file
    with nc.Dataset(STATIC_FILE) as ds:
        dom_soil = ds.variables['dominant_soiltype'][:] # (121, 181)
    dom_soil_flat = dom_soil.flatten()
    
    cols = [
        'pixel_id', 'date', 'month', 'season', 'pass', 'smap_sm', 'ssm_ol_daily',
        'irrigfrac_gmia', 'gmia_class', 'spatial_block_id',
        'qc_lis_rain', 'qc_lis_frozen', 'qc_lis_forest', 'qc_lis_snow',
        'qc_lis_high_veg', 'qc_lis_sat_limit', 'qc_lis_vwc_gt_5'
    ]
    
    dfs = []
    for y in [2016, 2017, 2018, 2019]:
        p = os.path.join(TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(p, columns=cols)
        # Filter observations passing non-wilt QC rules
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
    
    # Sort strictly by pixel_id, pass, day_idx
    full_df = full_df.sort_values(['pixel_id', 'pass', 'day_idx']).reset_index(drop=True)
    
    # Successive pairs check
    same_group = (full_df['pixel_id'] == full_df['pixel_id'].shift(1)) & (full_df['pass'] == full_df['pass'].shift(1))
    d_days = full_df['day_idx'] - full_df['day_idx'].shift(1)
    valid_pair = same_group & (d_days >= 1) & (d_days <= 3)
    
    t2 = full_df['day_idx'].values
    t1 = full_df['day_idx'].shift(1).fillna(-999).astype(int).values
    p_ids = full_df['pixel_id'].values
    dt = d_days.fillna(-999).astype(int).values
    
    # Precompute cumulative rain between t1 and t2 (including day t2 to ensure zero rain throughout)
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
    full_df['delta_ol'] = full_df['ssm_ol_daily'] - full_df['ssm_ol_daily'].shift(1)
    
    # Texture mapping
    textures = [TEXTURE_NAMES.get(dom_soil_flat[pix], "Other") for pix in full_df['pixel_id'].values]
    full_df['texture_name'] = textures
    
    # Subset to valid pairs only
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
    pairs_df['is_wet_ol'] = (pairs_df['delta_ol'] > sm_thresh) & (pairs_df['rain_cum_imerg'] < rain_thresh)
    
    rows = []
    np.random.seed(42)
    
    # Strata to evaluate:
    # 1. By GMIA class x Month (overall texture)
    # 2. By GMIA class x Season (overall texture)
    # 3. By GMIA class x Texture x JJA
    # 4. By GMIA class x Texture (all seasons)
    
    # Aggregate block-level stats
    gmia_classes = ['0%', '0-10%', '10-30%', '>30%', 'All']
    textures = ['All', 'Clay Loam', 'Loam', 'Sandy Clay Loam', 'Sandy Loam', 'Clay']
    months = list(range(1, 13)) + ['JJA', 'All']
    
    for gmia in gmia_classes:
        sub_g = pairs_df if gmia == 'All' else pairs_df[pairs_df['gmia_class'] == gmia]
        
        for tex in textures:
            sub_t = sub_g if tex == 'All' else sub_g[sub_g['texture_name'] == tex]
            
            for m in months:
                if m == 'All':
                    sub = sub_t
                elif m == 'JJA':
                    sub = sub_t[sub_t['season'] == 'JJA']
                else:
                    sub = sub_t[sub_t['month'] == m]
                    
                n_pairs = len(sub)
                if n_pairs == 0:
                    continue
                n_pixels = int(sub['pixel_id'].nunique())
                n_blocks = int(sub['spatial_block_id'].nunique())
                
                # Block aggregation
                block_stats = sub.groupby('spatial_block_id').agg(
                    ev_smap=('is_wet_smap', 'sum'),
                    ev_ol=('is_wet_ol', 'sum'),
                    n_tot=('valid_pair', 'count')
                )
                
                ev_smap_arr = block_stats['ev_smap'].values
                ev_ol_arr = block_stats['ev_ol'].values
                tot_arr = block_stats['n_tot'].values
                
                smap_rate, smap_ci_l, smap_ci_h = compute_stratum_bootstrap_rate(ev_smap_arr, tot_arr)
                ol_rate, ol_ci_l, ol_ci_h = compute_stratum_bootstrap_rate(ev_ol_arr, tot_arr)
                
                rows.append({
                    'gmia_class': gmia,
                    'texture': tex,
                    'temporal_scope': str(m),
                    'n_pairs': n_pairs,
                    'n_pixels': n_pixels,
                    'n_blocks': n_blocks,
                    'smap_wet_events': int(ev_smap_arr.sum()),
                    'smap_freq_pct': smap_rate,
                    'smap_ci_low': smap_ci_l,
                    'smap_ci_high': smap_ci_h,
                    'ol_wet_events': int(ev_ol_arr.sum()),
                    'ol_freq_pct': ol_rate,
                    'ol_ci_low': ol_ci_l,
                    'ol_ci_high': ol_ci_h,
                    'smap_ol_diff_pct': smap_rate - ol_rate
                })
                
    return pd.DataFrame(rows)


def plot_figure_f6(df_diag):
    """Generate Figure F6 (FR and EN) showing wetting frequencies."""
    print("Generating Figure F6 (Dynamic Irrigation Wetting Events)...")
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(FIG_EN_DIR, exist_ok=True)
    
    # Prepare monthly data for Clay Loam
    cl_months = df_diag[(df_diag['texture'] == 'Clay Loam') & 
                        (df_diag['temporal_scope'].isin([str(i) for i in range(1, 13)]))].copy()
    cl_months['month_int'] = cl_months['temporal_scope'].astype(int)
    cl_months = cl_months.sort_values('month_int')
    
    # Prepare JJA data across GMIA classes for Clay Loam and Loam
    jja_data = df_diag[(df_diag['temporal_scope'] == 'JJA') & 
                       (df_diag['gmia_class'] != 'All') & 
                       (df_diag['texture'].isin(['Clay Loam', 'Loam']))].copy()
                       
    gmia_order = ['0%', '0-10%', '10-30%', '>30%']
    
    for lang, out_dir, fname in [
        ('fr', FIG_DIR, 'F6_dynamic_irrigation_wetting_events.png'),
        ('en', FIG_EN_DIR, 'F6_dynamic_irrigation_wetting_events_en.png')
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
        
        # Panel A: Monthly cycle for Clay Loam (GMIA >30% vs GMIA 0%)
        ax0 = axes[0]
        months_arr = np.arange(1, 13)
        
        # SMAP >30%
        sub_gt30 = cl_months[cl_months['gmia_class'] == '>30%'].set_index('month_int').reindex(months_arr)
        # SMAP 0%
        sub_0 = cl_months[cl_months['gmia_class'] == '0%'].set_index('month_int').reindex(months_arr)
        # OL control >30%
        ol_gt30 = cl_months[cl_months['gmia_class'] == '>30%'].set_index('month_int').reindex(months_arr)
        
        lbl_smap_gt30 = "SMAP (GMIA >30%)" if lang == 'fr' else "SMAP (GMIA >30%)"
        lbl_smap_0 = "SMAP (GMIA 0%)" if lang == 'fr' else "SMAP (GMIA 0%)"
        lbl_ol = "Noah-MP OL témoin (GMIA >30%)" if lang == 'fr' else "Noah-MP OL control (GMIA >30%)"
        
        ax0.plot(months_arr, sub_gt30['smap_freq_pct'].values, 'o-', color='#1b9e77', linewidth=2.2, label=lbl_smap_gt30)
        ax0.fill_between(months_arr, sub_gt30['smap_ci_low'].values, sub_gt30['smap_ci_high'].values, color='#1b9e77', alpha=0.2)
        
        ax0.plot(months_arr, sub_0['smap_freq_pct'].values, 's--', color='#7570b3', linewidth=1.8, label=lbl_smap_0)
        ax0.fill_between(months_arr, sub_0['smap_ci_low'].values, sub_0['smap_ci_high'].values, color='#7570b3', alpha=0.15)
        
        ax0.plot(months_arr, ol_gt30['ol_freq_pct'].values, '^:', color='#d95f02', linewidth=1.8, label=lbl_ol)
        
        ax0.set_xlabel("Mois" if lang == 'fr' else "Month", fontsize=11, fontweight='bold')
        ax0.set_ylabel("Fréquence d'humidification sans pluie (%)" if lang == 'fr' else "Wetting frequency without rain (%)", fontsize=11, fontweight='bold')
        ax0.set_title("(a) Cycle mensuel en Clay Loam (ΔSM > 0.03 m³/m³, Pluie < 1 mm)" if lang == 'fr' else "(a) Monthly cycle in Clay Loam (ΔSM > 0.03 m³/m³, Rain < 1 mm)", fontsize=11)
        ax0.set_xticks(months_arr)
        ax0.grid(True, linestyle='--', alpha=0.5)
        ax0.legend(frameon=True, fontsize=9.5)
        
        # Panel B: JJA frequency by GMIA class
        ax1 = axes[1]
        x_indices = np.arange(len(gmia_order))
        width = 0.35
        
        # Clay Loam JJA
        cl_jja = jja_data[jja_data['texture'] == 'Clay Loam'].set_index('gmia_class').reindex(gmia_order)
        # Loam JJA
        loam_jja = jja_data[jja_data['texture'] == 'Loam'].set_index('gmia_class').reindex(gmia_order)
        
        yerr_cl = [(cl_jja['smap_freq_pct'] - cl_jja['smap_ci_low']).values, (cl_jja['smap_ci_high'] - cl_jja['smap_freq_pct']).values]
        yerr_loam = [(loam_jja['smap_freq_pct'] - loam_jja['smap_ci_low']).values, (loam_jja['smap_ci_high'] - loam_jja['smap_freq_pct']).values]
        
        ax1.bar(x_indices - width/2, cl_jja['smap_freq_pct'].values, width, yerr=yerr_cl, capsize=4, color='#2b83ba', label='Clay Loam (SMAP)')
        ax1.bar(x_indices + width/2, loam_jja['smap_freq_pct'].values, width, yerr=yerr_loam, capsize=4, color='#fdae61', label='Loam (SMAP)')
        
        # Add OL points on top
        ax1.plot(x_indices - width/2, cl_jja['ol_freq_pct'].values, 'ro', markersize=6, label='Clay Loam (OL)')
        ax1.plot(x_indices + width/2, loam_jja['ol_freq_pct'].values, 'mo', markersize=6, label='Loam (OL)')
        
        ax1.set_xlabel("Classe d'irrigation GMIA" if lang == 'fr' else "GMIA Irrigation Class", fontsize=11, fontweight='bold')
        ax1.set_ylabel("Fréquence d'humidification JJA (%)" if lang == 'fr' else "JJA Wetting frequency (%)", fontsize=11, fontweight='bold')
        ax1.set_title("(b) Signal d'été (JJA) par classe GMIA" if lang == 'fr' else "(b) Summer signal (JJA) by GMIA class", fontsize=11)
        ax1.set_xticks(x_indices)
        ax1.set_xticklabels(gmia_order, fontsize=10)
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax1.legend(frameon=True, fontsize=9.5)
        
        plt.tight_layout()
        fpath = os.path.join(out_dir, fname)
        plt.savefig(fpath, dpi=300)
        plt.close()
        print(f"Saved: {fpath}")


def generate_markdown_report(df_default, df_sens, pairs_df, thresh_cfg):
    """Generate Phase 2a-bis Markdown report with zero hardcoded numbers and zero interpretive words."""
    print("Generating Phase 2a-bis Markdown report...")
    out_md = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic_report.md")
    
    n_pairs_tot = len(pairs_df)
    n_pixels_tot = int(pairs_df['pixel_id'].nunique())
    n_blocks_tot = int(pairs_df['spatial_block_id'].nunique())
    
    # Overall numbers
    ov_default = df_default[(df_default['gmia_class'] == 'All') & 
                            (df_default['texture'] == 'All') & 
                            (df_default['temporal_scope'] == 'All')].iloc[0]
                            
    # JJA numbers for GMIA classes in Clay Loam
    cl_jja_rows = []
    for gmia in ['0%', '0-10%', '10-30%', '>30%']:
        sub = df_default[(df_default['texture'] == 'Clay Loam') & 
                         (df_default['temporal_scope'] == 'JJA') & 
                         (df_default['gmia_class'] == gmia)]
        if len(sub) > 0:
            r = sub.iloc[0]
            cl_jja_rows.append(
                f"| {r['gmia_class']} | {r['n_pairs']:,} | {r['n_pixels']:,} | {r['n_blocks']:,} | "
                f"{r['smap_wet_events']:,} | {r['smap_freq_pct']:.3f}% [{r['smap_ci_low']:.3f}%, {r['smap_ci_high']:.3f}%] | "
                f"{r['ol_wet_events']:,} | {r['ol_freq_pct']:.3f}% [{r['ol_ci_low']:.3f}%, {r['ol_ci_high']:.3f}%] | "
                f"{r['smap_ol_diff_pct']:+.3f}% |"
            )
    table_cl_jja = "\n".join(cl_jja_rows)
    
    # Loam JJA rows
    loam_jja_rows = []
    for gmia in ['0%', '0-10%', '10-30%', '>30%']:
        sub = df_default[(df_default['texture'] == 'Loam') & 
                         (df_default['temporal_scope'] == 'JJA') & 
                         (df_default['gmia_class'] == gmia)]
        if len(sub) > 0:
            r = sub.iloc[0]
            loam_jja_rows.append(
                f"| {r['gmia_class']} | {r['n_pairs']:,} | {r['n_pixels']:,} | {r['n_blocks']:,} | "
                f"{r['smap_wet_events']:,} | {r['smap_freq_pct']:.3f}% [{r['smap_ci_low']:.3f}%, {r['smap_ci_high']:.3f}%] | "
                f"{r['ol_wet_events']:,} | {r['ol_freq_pct']:.3f}% [{r['ol_ci_low']:.3f}%, {r['ol_ci_high']:.3f}%] | "
                f"{r['smap_ol_diff_pct']:+.3f}% |"
            )
    table_loam_jja = "\n".join(loam_jja_rows)
    
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
            r_gt = sub_gt30.iloc[0]
            r_0 = sub_0.iloc[0]
            month_rows.append(
                f"| {m:02d} | {r_gt['smap_freq_pct']:.3f}% [{r_gt['smap_ci_low']:.3f}%, {r_gt['smap_ci_high']:.3f}%] | "
                f"{r_0['smap_freq_pct']:.3f}% [{r_0['smap_ci_low']:.3f}%, {r_0['smap_ci_high']:.3f}%] | "
                f"{r_gt['smap_freq_pct'] - r_0['smap_freq_pct']:+.3f}% | "
                f"{r_gt['ol_freq_pct']:.3f}% | {r_0['ol_freq_pct']:.3f}% | "
                f"{r_gt['n_pixels']:,} / {r_0['n_pixels']:,} | {r_gt['n_blocks']:,} / {r_0['n_blocks']:,} |"
            )
    table_months = "\n".join(month_rows)
    
    # Sensitivity summary table
    sens_rows = []
    for _, r_s in df_sens.iterrows():
        sens_rows.append(
            f"| {r_s['threshold_config']} | ΔSM > {r_s['delta_sm_min']:.2f} m³/m³, Rain < {r_s['rain_imerg_max']:.1f} mm | "
            f"{r_s['jja_gt30_smap_pct']:.3f}% | {r_s['jja_0_smap_pct']:.3f}% | "
            f"{r_s['jja_diff_smap_pct']:+.3f}% | {r_s['jja_gt30_ol_pct']:.3f}% |"
        )
    table_sens = "\n".join(sens_rows)
    
    report_md = f"""# Rapport de Diagnostic Phase 2a-bis : Signal d'Irrigation Dynamique

**Période d'évaluation** : 2016-01-01 au 2019-12-31 (4 années complètes).  
**Statut de l'année 2020** : 2020 reste strictement scellée.  
**Population d'observation** : Paires d'observations SMAP successives du même passage orbital (AM avec AM, PM avec PM), espacées de 1 à 3 jours, issues des observations conformes aux règles QC hors flétrissement.  
**Effectif total analysé** : {n_pairs_tot:,} paires successives ({n_pixels_tot:,} mailles distinctes, {n_blocks_tot:,} blocs spatiaux distincts de 1°).

---

## 1. Définition Opérationnelle et Seuils

Un événement d'humidification sans précipitation est défini par :
- Un intervalle temporel de 1 à 3 jours entre deux observations du même passage au même pixel.
- Une élévation d'humidité du sol : $\Delta \text{{SM}} > \tau_{{\text{{sm}}}}$.
- Un cumul de précipitations IMERG sur l'intervalle : $\text{{Rain}}_{{t_1 \to t_2}} < \tau_{{\text{{rain}}}}$.
- Configuration par défaut : $\tau_{{\text{{sm}}}} = {thresh_cfg['default']['delta_sm_min']:.2f}\text{{ m}}^3/\text{{m}}^3$, $\tau_{{\text{{rain}}}} = {thresh_cfg['default']['rain_imerg_max']:.1f}\text{{ mm}}$.
- Le modèle Open Loop (Noah-MP sans irrigation) subit exactement le même calcul sur les mêmes paires d'observations pour servir de témoin négatif.

---

## 2. Statistiques Globales (Configuration par Défaut)

- **Nombre total de paires de passages successifs** : {n_pairs_tot:,}
- **Événements détectés sur SMAP (ensemble du domaine et de l'année)** : {ov_default['smap_wet_events']:,} ({ov_default['smap_freq_pct']:.3f}%, IC 95% : [{ov_default['smap_ci_low']:.3f}%, {ov_default['smap_ci_high']:.3f}%])
- **Événements détectés sur l'Open Loop (témoin)** : {ov_default['ol_wet_events']:,} ({ov_default['ol_freq_pct']:.3f}%, IC 95% : [{ov_default['ol_ci_low']:.3f}%, {ov_default['ol_ci_high']:.3f}%])
- **Écart global SMAP − OL** : {ov_default['smap_ol_diff_pct']:+.3f}%

---

## 3. Signal d'Été (JJA) par Classe d'Irrigation GMIA à Texture Fixée

### 3.1. Texture Clay Loam (JJA)

| Classe GMIA | Paires | Pixels | Blocs 1° | Événements SMAP | Fréquence SMAP [IC 95%] | Événements OL | Fréquence OL [IC 95%] | Contraste SMAP − OL |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{table_cl_jja}

### 3.2. Texture Loam (JJA)

| Classe GMIA | Paires | Pixels | Blocs 1° | Événements SMAP | Fréquence SMAP [IC 95%] | Événements OL | Fréquence OL [IC 95%] | Contraste SMAP − OL |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{table_loam_jja}

---

## 4. Cycle Mensuel en Clay Loam : Contraste GMIA >30 % vs 0 %

| Mois | SMAP >30% [IC 95%] | SMAP 0% [IC 95%] | Contraste SMAP (>30% − 0%) | OL >30% | OL 0% | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
{table_months}

---

## 5. Analyse de Sensibilité aux Seuils

| Configuration | Seuils appliqués | SMAP JJA >30% | SMAP JJA 0% | Contraste SMAP JJA | OL JJA >30% |
|:---|:---|:---:|:---:|:---:|:---:|
{table_sens}

---

## 6. Référence aux Figures et Données Produites

- **Figure F6** : `figures/descriptive/F6_dynamic_irrigation_wetting_events.png` (version EN : `figures/descriptive/en/F6_dynamic_irrigation_wetting_events_en.png`)
- **Table complète par strate** : `reports/phase2a_irrigation_diagnostic.csv`
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
    lax_cfg = thresh_cfg['sensitivity_lax']
    strict_cfg = thresh_cfg['sensitivity_strict']
    
    imerg_flat = load_continuous_imerg()
    pairs_df = load_and_pair_observations(imerg_flat)
    
    # 1. Default evaluation
    print(f"Evaluating default thresholds: delta_sm > {def_cfg['delta_sm_min']}, rain < {def_cfg['rain_imerg_max']} mm...")
    df_default = evaluate_diagnostic(pairs_df, def_cfg['delta_sm_min'], def_cfg['rain_imerg_max'])
    csv_default = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic.csv")
    df_default.to_csv(csv_default, index=False)
    print(f"Saved default diagnostic CSV: {csv_default}")
    
    # 2. Sensitivity evaluations
    sens_rows = []
    for cfg_key in ['default', 'sensitivity_lax', 'sensitivity_strict']:
        cfg = thresh_cfg[cfg_key]
        if cfg_key == 'default':
            df_s = df_default
        else:
            print(f"Evaluating sensitivity configuration: {cfg_key}...")
            df_s = evaluate_diagnostic(pairs_df, cfg['delta_sm_min'], cfg['rain_imerg_max'])
        
        # Extract JJA Clay Loam >30% and 0%
        sub_gt = df_s[(df_s['texture'] == 'Clay Loam') & (df_s['temporal_scope'] == 'JJA') & (df_s['gmia_class'] == '>30%')]
        sub_0 = df_s[(df_s['texture'] == 'Clay Loam') & (df_s['temporal_scope'] == 'JJA') & (df_s['gmia_class'] == '0%')]
        
        pct_gt = float(sub_gt['smap_freq_pct'].iloc[0]) if len(sub_gt) > 0 else np.nan
        pct_0 = float(sub_0['smap_freq_pct'].iloc[0]) if len(sub_0) > 0 else np.nan
        ol_gt = float(sub_gt['ol_freq_pct'].iloc[0]) if len(sub_gt) > 0 else np.nan
        ol_0 = float(sub_0['ol_freq_pct'].iloc[0]) if len(sub_0) > 0 else np.nan
        
        sens_rows.append({
            'threshold_config': cfg_key,
            'delta_sm_min': cfg['delta_sm_min'],
            'rain_imerg_max': cfg['rain_imerg_max'],
            'jja_gt30_smap_pct': pct_gt,
            'jja_0_smap_pct': pct_0,
            'jja_diff_smap_pct': pct_gt - pct_0,
            'jja_gt30_ol_pct': ol_gt,
            'jja_0_ol_pct': ol_0
        })
        
    df_sens = pd.DataFrame(sens_rows)
    csv_sens = os.path.join(REPORTS_DIR, "phase2a_irrigation_diagnostic_sensitivity.csv")
    df_sens.to_csv(csv_sens, index=False)
    print(f"Saved sensitivity CSV: {csv_sens}")
    
    # 3. Generate Figure F6
    plot_figure_f6(df_default)
    
    # 4. Generate Markdown report
    generate_markdown_report(df_default, df_sens, pairs_df, thresh_cfg)
    
    print(f"=== PHASE 2a-bis COMPLETED IN {time.time()-t_start:.1f}s ===")


if __name__ == '__main__':
    main()
