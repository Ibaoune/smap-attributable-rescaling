#!/usr/bin/env python3
"""
06_generate_descriptive_figures.py

Generate descriptive figures F1 to F4 and their associated source CSV tables:
  F1. Observation frequency per pixel (AM and PM, before/after QC filter).
  F2. Seasonal bias map (d = SMAP - SSM_OL) for AM and PM.
  F3. d as a function of GMIA irrigation class, stratified by soil texture and season,
      with 95% confidence intervals via spatial block bootstrap (1° blocks) and distinct pixel counts.
  F4. Share of observations rejected by wilting threshold: spatial map and monthly cycle.

Generates both:
  - French versions in figures/descriptive/
  - English versions in figures/descriptive/en/ (purely descriptive, zero interpretation phrases)
"""

import os
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT_TABLE_DIR = os.path.join(BASE_DIR, "data", "processed", "modeling_table")
STATIC_FILE = os.path.join(BASE_DIR, "data", "interim", "static", "static_predictors.nc")
OUT_FIG_DIR = os.path.join(BASE_DIR, "figures", "descriptive")
OUT_FIG_DIR_EN = os.path.join(BASE_DIR, "figures", "descriptive", "en")

SOIL_NAMES = {
    1: "Sand", 2: "Loamy Sand", 3: "Sandy Loam", 4: "Silt Loam",
    5: "Silt", 6: "Loam", 7: "Sandy Clay Loam", 8: "Silty Clay Loam",
    9: "Clay Loam", 10: "Sandy Clay", 11: "Silty Clay", 12: "Clay",
    13: "Organic Material", 14: "Water", 15: "Bedrock", 16: "Other"
}


def load_all_data():
    """Load static predictors and all 5 annual modeling tables."""
    print("Loading static predictors...")
    with netCDF4.Dataset(STATIC_FILE) as ds:
        lat_1d = ds.variables['lat'][:]
        lon_1d = ds.variables['lon'][:]
        landmask = ds.variables['landmask'][:]
        dom_soil = ds.variables['dominant_soiltype'][:]

    print("Loading annual modeling tables (2016-2020)...")
    dfs = []
    years = [2016, 2017, 2018, 2019, 2020]
    for y in years:
        pq = os.path.join(OUT_TABLE_DIR, f"sar_modeling_table_{y}.parquet")
        df_y = pd.read_parquet(pq, columns=[
            'year', 'month', 'pass', 'pixel_id', 'spatial_block_id',
            'season', 'gmia_class', 'target_d', 'qc_lis_combined', 'qc_lis_wilt_limit'
        ])
        print(f"  Loaded {pq}: {len(df_y):,} rows")
        dfs.append(df_y)
    df_all = pd.concat(dfs, ignore_index=True)
    print(f"Total observations loaded: {len(df_all):,}")
    return lat_1d, lon_1d, landmask, dom_soil, df_all


def generate_f1(lat_1d, lon_1d, landmask, df_all):
    """F1: Observation frequency per pixel, AM and PM, before and after LIS masks."""
    print("--- Generating Figure F1 & Source CSV ---")
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    
    raw_am_grid = np.zeros((n_lat, n_lon), dtype=np.int32)
    filt_am_grid = np.zeros((n_lat, n_lon), dtype=np.int32)
    raw_pm_grid = np.zeros((n_lat, n_lon), dtype=np.int32)
    filt_pm_grid = np.zeros((n_lat, n_lon), dtype=np.int32)
    
    df_am = df_all[df_all['pass'] == 'AM']
    df_pm = df_all[df_all['pass'] == 'PM']
    
    am_raw_cnt = df_am.groupby('pixel_id').size()
    am_filt_cnt = df_am[~df_am['qc_lis_combined']].groupby('pixel_id').size()
    pm_raw_cnt = df_pm.groupby('pixel_id').size()
    pm_filt_cnt = df_pm[~df_pm['qc_lis_combined']].groupby('pixel_id').size()
    
    for pid, c in am_raw_cnt.items():
        raw_am_grid[pid // n_lon, pid % n_lon] = c
    for pid, c in am_filt_cnt.items():
        filt_am_grid[pid // n_lon, pid % n_lon] = c
    for pid, c in pm_raw_cnt.items():
        raw_pm_grid[pid // n_lon, pid % n_lon] = c
    for pid, c in pm_filt_cnt.items():
        filt_pm_grid[pid // n_lon, pid % n_lon] = c
        
    mask_sea = (landmask == 0)
    raw_am_m = np.ma.masked_where(mask_sea, raw_am_grid)
    filt_am_m = np.ma.masked_where(mask_sea, filt_am_grid)
    raw_pm_m = np.ma.masked_where(mask_sea, raw_pm_grid)
    filt_pm_m = np.ma.masked_where(mask_sea, filt_pm_grid)
    
    # Save CSV
    csv_rows = []
    for r in range(n_lat):
        for c in range(n_lon):
            if landmask[r, c] == 1:
                pid = r * n_lon + c
                r_am = int(raw_am_grid[r, c])
                f_am = int(filt_am_grid[r, c])
                r_pm = int(raw_pm_grid[r, c])
                f_pm = int(filt_pm_grid[r, c])
                csv_rows.append({
                    'pixel_id': pid,
                    'lat': lat_1d[r],
                    'lon': lon_1d[c],
                    'count_raw_am': r_am,
                    'count_filtered_am': f_am,
                    'pct_retained_am': (f_am / r_am * 100.0) if r_am > 0 else 0.0,
                    'count_raw_pm': r_pm,
                    'count_filtered_pm': f_pm,
                    'pct_retained_pm': (f_pm / r_pm * 100.0) if r_pm > 0 else 0.0,
                })
    df_csv = pd.DataFrame(csv_rows)
    csv_path = os.path.join(OUT_FIG_DIR, "F1_observation_frequency.csv")
    df_csv.to_csv(csv_path, index=False)
    print(f"  Saved F1 source CSV: {csv_path}")
    
    extent = [lon_1d.min(), lon_1d.max(), lat_1d.min(), lat_1d.max()]
    vmax = max(raw_am_grid.max(), raw_pm_grid.max())
    
    # French version
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    panels_fr = [
        (axes[0, 0], raw_am_m, "AM - Avant masques LIS (Brut)", "Blues"),
        (axes[0, 1], filt_am_m, "AM - Après masques LIS (Filtré)", "Blues"),
        (axes[1, 0], raw_pm_m, "PM - Avant masques LIS (Brut)", "Purples"),
        (axes[1, 1], filt_pm_m, "PM - Après masques LIS (Filtré)", "Purples"),
    ]
    for ax, data, title, cmap in panels_fr:
        im = ax.imshow(data, origin='lower', extent=extent, cmap=cmap, vmin=0, vmax=vmax)
        ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
        ax.set_xlabel("Longitude (°E)", fontsize=10)
        ax.set_ylabel("Latitude (°N)", fontsize=10)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Nombre d'observations (2016-2020)", fontsize=9)
        ax.grid(True, linestyle='--', alpha=0.4)
    plt.suptitle("F1. Fréquence spatiale des observations SMAP (2016-2020) : Passages AM/PM avant et après masques LIS",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig_path_fr = os.path.join(OUT_FIG_DIR, "F1_observation_frequency.png")
    plt.savefig(fig_path_fr, dpi=300)
    plt.close()
    
    # English version
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    panels_en = [
        (axes[0, 0], raw_am_m, "AM - Raw (Before LIS QC)", "Blues"),
        (axes[0, 1], filt_am_m, "AM - Filtered (After LIS QC)", "Blues"),
        (axes[1, 0], raw_pm_m, "PM - Raw (Before LIS QC)", "Purples"),
        (axes[1, 1], filt_pm_m, "PM - Filtered (After LIS QC)", "Purples"),
    ]
    for ax, data, title, cmap in panels_en:
        im = ax.imshow(data, origin='lower', extent=extent, cmap=cmap, vmin=0, vmax=vmax)
        ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
        ax.set_xlabel("Longitude (°E)", fontsize=10)
        ax.set_ylabel("Latitude (°N)", fontsize=10)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Observation count (2016-2020)", fontsize=9)
        ax.grid(True, linestyle='--', alpha=0.4)
    plt.suptitle("F1. Spatial observation count of SMAP (2016-2020): AM and PM passes before and after LIS QC masks",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig_path_en = os.path.join(OUT_FIG_DIR_EN, "F1_observation_frequency_en.png")
    plt.savefig(fig_path_en, dpi=300)
    plt.close()
    print(f"  Saved Figure F1 (FR & EN): {fig_path_fr}, {fig_path_en}")


def generate_f2(lat_1d, lon_1d, landmask, df_all):
    """F2: Map of seasonal mean d (AM and PM separated)."""
    print("--- Generating Figure F2 & Source CSV ---")
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    passes = ['AM', 'PM']
    
    grids = {f"{p}_{s}": np.full((n_lat, n_lon), np.nan, dtype=np.float32)
             for p in passes for s in seasons}
    
    grouped = df_all.groupby(['pass', 'season', 'pixel_id'])['target_d'].mean()
    for (p, s, pid), val in grouped.items():
        grids[f"{p}_{s}"][pid // n_lon, pid % n_lon] = val
        
    csv_rows = []
    for r in range(n_lat):
        for c in range(n_lon):
            if landmask[r, c] == 1:
                pid = r * n_lon + c
                row = {'pixel_id': pid, 'lat': lat_1d[r], 'lon': lon_1d[c]}
                for p in passes:
                    for s in seasons:
                        v = grids[f"{p}_{s}"][r, c]
                        row[f"mean_d_{p.lower()}_{s.lower()}"] = np.round(v, 5) if not np.isnan(v) else np.nan
                csv_rows.append(row)
    df_csv = pd.DataFrame(csv_rows)
    csv_path = os.path.join(OUT_FIG_DIR, "F2_season_bias_maps.csv")
    df_csv.to_csv(csv_path, index=False)
    print(f"  Saved F2 source CSV: {csv_path}")
    
    extent = [lon_1d.min(), lon_1d.max(), lat_1d.min(), lat_1d.max()]
    mask_sea = (landmask == 0)
    vlim = 0.12
    
    # French version
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True, sharey=True)
    for r_idx, p in enumerate(passes):
        for c_idx, s in enumerate(seasons):
            ax = axes[r_idx, c_idx]
            g = np.ma.masked_where(mask_sea, grids[f"{p}_{s}"])
            im = ax.imshow(g, origin='lower', extent=extent, cmap='coolwarm', vmin=-vlim, vmax=vlim)
            ax.set_title(f"{p} - {s}", fontsize=12, fontweight='bold', pad=8)
            if r_idx == 1:
                ax.set_xlabel("Longitude (°E)", fontsize=10)
            if c_idx == 0:
                ax.set_ylabel(f"Latitude (°N)\n[{p}]", fontsize=11, fontweight='bold')
            ax.grid(True, linestyle='--', alpha=0.4)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.70, 0.025])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Moyenne de $d = \\mathrm{SMAP} - \\mathrm{SSM}_{\\mathrm{OL}}$ [$\\mathrm{m}^3/\\mathrm{m}^3$]", fontsize=11)
    plt.suptitle("F2. Cartographie saisonnière de l'écart moyen $d$ (2016-2020) : Passages AM (06h) et PM (18h)",
                 fontsize=15, fontweight='bold', y=0.98)
    plt.subplots_adjust(top=0.92, bottom=0.14, hspace=0.18, wspace=0.12)
    fig_path_fr = os.path.join(OUT_FIG_DIR, "F2_season_bias_maps.png")
    plt.savefig(fig_path_fr, dpi=300)
    plt.close()

    # English version
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True, sharey=True)
    for r_idx, p in enumerate(passes):
        for c_idx, s in enumerate(seasons):
            ax = axes[r_idx, c_idx]
            g = np.ma.masked_where(mask_sea, grids[f"{p}_{s}"])
            im = ax.imshow(g, origin='lower', extent=extent, cmap='coolwarm', vmin=-vlim, vmax=vlim)
            ax.set_title(f"{p} - {s}", fontsize=12, fontweight='bold', pad=8)
            if r_idx == 1:
                ax.set_xlabel("Longitude (°E)", fontsize=10)
            if c_idx == 0:
                ax.set_ylabel(f"Latitude (°N)\n[{p}]", fontsize=11, fontweight='bold')
            ax.grid(True, linestyle='--', alpha=0.4)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.70, 0.025])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Mean difference $d = \\mathrm{SMAP} - \\mathrm{SSM}_{\\mathrm{OL}}$ [$\\mathrm{m}^3/\\mathrm{m}^3$]", fontsize=11)
    plt.suptitle("F2. Seasonal maps of mean difference d = SMAP - SSM_OL (2016-2020): AM (06:00) and PM (18:00) passes",
                 fontsize=15, fontweight='bold', y=0.98)
    plt.subplots_adjust(top=0.92, bottom=0.14, hspace=0.18, wspace=0.12)
    fig_path_en = os.path.join(OUT_FIG_DIR_EN, "F2_season_bias_maps_en.png")
    plt.savefig(fig_path_en, dpi=300)
    plt.close()
    print(f"  Saved Figure F2 (FR & EN): {fig_path_fr}, {fig_path_en}")


def generate_f3(dom_soil, df_all):
    """F3: d vs GMIA class stratified by dominant soil texture and season (JJA separate),
    with 95% CIs via spatial block bootstrap on 1° blocks, and distinct pixel counts for each point.
    """
    print("--- Generating Figure F3 & Source CSV ---")
    n_lon = 181
    
    pids = df_all['pixel_id'].values
    r_idx = pids // n_lon
    c_idx = pids % n_lon
    obs_soil = dom_soil[r_idx, c_idx]
    df_all['dominant_soiltype'] = obs_soil
    
    top_soils = [6, 9, 3, 7, 12]  # Loam, Clay Loam, Sandy Loam, Sandy Clay Loam, Clay
    df_sub = df_all[df_all['dominant_soiltype'].isin(top_soils)].copy()
    df_sub['soil_name'] = df_sub['dominant_soiltype'].map(SOIL_NAMES)
    
    df_sub['season_group_fr'] = np.where(df_sub['season'] == 'JJA', 'JJA (Été)', 'Autres saisons (DJF, MAM, SON)')
    df_sub['season_group_en'] = np.where(df_sub['season'] == 'JJA', 'JJA (Summer)', 'Other seasons (DJF, MAM, SON)')
    
    gmia_order = ['0%', '0-10%', '10-30%', '>30%']
    
    # Compute stratified statistics with 1° spatial block bootstrap CIs & distinct pixel counts
    rng = np.random.RandomState(42)
    records = []
    
    for (sg_fr, sg_en, soil, g_cls), group_df in df_sub.groupby(['season_group_fr', 'season_group_en', 'soil_name', 'gmia_class']):
        n_obs = len(group_df)
        n_pix = int(group_df['pixel_id'].nunique())
        m_val = float(group_df['target_d'].mean())
        s_val = float(group_df['target_d'].std())
        
        # 1° spatial block bootstrap (B=500)
        blocks = group_df['spatial_block_id'].unique()
        if len(blocks) > 1:
            block_grp = group_df.groupby('spatial_block_id')['target_d'].agg(['mean', 'count'])
            b_means = block_grp['mean'].values
            b_counts = block_grp['count'].values
            n_b = len(b_means)
            boot_means = np.empty(500, dtype=np.float64)
            for b in range(500):
                sampled_idx = rng.choice(n_b, size=n_b, replace=True)
                w = b_counts[sampled_idx]
                boot_means[b] = np.sum(w * b_means[sampled_idx]) / np.sum(w)
            ci_lower = float(np.percentile(boot_means, 2.5))
            ci_upper = float(np.percentile(boot_means, 97.5))
            ci_half = (ci_upper - ci_lower) / 2.0
        else:
            ci_half = 1.96 * s_val / np.sqrt(max(1, n_pix))
            ci_lower = m_val - ci_half
            ci_upper = m_val + ci_half
            
        records.append({
            'season_group_fr': sg_fr,
            'season_group_en': sg_en,
            'soil_name': soil,
            'gmia_class': g_cls,
            'mean': m_val,
            'std': s_val,
            'n_obs': n_obs,
            'n_distinct_pixels': n_pix,
            'ci95_half': ci_half,
            'ci95_lower': ci_lower,
            'ci95_upper': ci_upper
        })
        
    stats = pd.DataFrame(records)
    csv_path = os.path.join(OUT_FIG_DIR, "F3_d_by_gmia_texture_season.csv")
    stats.to_csv(csv_path, index=False)
    print(f"  Saved F3 source CSV with distinct pixel counts and block bootstrap CIs: {csv_path}")
    
    colors = {
        'Loam': '#1f77b4',
        'Clay Loam': '#d62728',
        'Sandy Loam': '#2ca02c',
        'Sandy Clay Loam': '#ff7f0e',
        'Clay': '#9467bd'
    }
    
    # 1. French Plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True)
    groups_fr = ['Autres saisons (DJF, MAM, SON)', 'JJA (Été)']
    panel_titles_fr = [
        "A. Saisons humides et intermédiaires (DJF, MAM, SON)",
        "B. Saison estivale (JJA, panneau dédié)"
    ]
    for ax_idx, (grp, title) in enumerate(zip(groups_fr, panel_titles_fr)):
        ax = axes[ax_idx]
        sub_grp = stats[stats['season_group_fr'] == grp]
        for soil in ['Loam', 'Clay Loam', 'Sandy Loam', 'Sandy Clay Loam', 'Clay']:
            s_data = sub_grp[sub_grp['soil_name'] == soil]
            if len(s_data) == 0:
                continue
            s_data = s_data.set_index('gmia_class').reindex(gmia_order).reset_index()
            x = np.arange(len(gmia_order))
            y = s_data['mean'].values
            ci = s_data['ci95_half'].values
            npix = s_data['n_distinct_pixels'].values
            
            ax.errorbar(x, y, yerr=ci, label=soil, fmt='-o', color=colors[soil],
                        capsize=5, capthick=1.5, linewidth=2, markersize=6)
            for xi, yi, pxi in zip(x, y, npix):
                if not np.isnan(yi):
                    ax.annotate(f"{pxi}px", (xi, yi), textcoords="offset points",
                                xytext=(0, 6 if soil in ['Loam', 'Clay Loam'] else -10),
                                ha='center', fontsize=6.5, color=colors[soil], alpha=0.85)

        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xticks(range(len(gmia_order)))
        ax.set_xticklabels(gmia_order, fontsize=10)
        ax.set_xlabel("Classe d'équipement d'irrigation GMIA (%)", fontsize=11)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.6, linewidth=1)
        ax.grid(True, linestyle=':', alpha=0.5)
        
    axes[0].set_ylabel("Écart moyen $d = \\mathrm{SMAP} - \\mathrm{SSM}_{\\mathrm{OL}}$ [$\\mathrm{m}^3/\\mathrm{m}^3$]", fontsize=11)
    axes[0].legend(title="Texture dominante", fontsize=9, title_fontsize=10, loc='lower left')
    plt.suptitle("F3. Écart $d$ en fonction de la classe GMIA, stratifié par texture dominante et par saison (IC 95% bootstrap 1° & nb pixels distincts)",
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig_path_fr = os.path.join(OUT_FIG_DIR, "F3_d_by_gmia_texture_season.png")
    plt.savefig(fig_path_fr, dpi=300)
    plt.close()
    
    # 2. English Plot (purely descriptive, zero interpretation phrases)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True)
    groups_en = ['Other seasons (DJF, MAM, SON)', 'JJA (Summer)']
    panel_titles_en = [
        "A. Wet and intermediate seasons (DJF, MAM, SON)",
        "B. Summer season (JJA)"
    ]
    for ax_idx, (grp, title) in enumerate(zip(groups_en, panel_titles_en)):
        ax = axes[ax_idx]
        sub_grp = stats[stats['season_group_en'] == grp]
        for soil in ['Loam', 'Clay Loam', 'Sandy Loam', 'Sandy Clay Loam', 'Clay']:
            s_data = sub_grp[sub_grp['soil_name'] == soil]
            if len(s_data) == 0:
                continue
            s_data = s_data.set_index('gmia_class').reindex(gmia_order).reset_index()
            x = np.arange(len(gmia_order))
            y = s_data['mean'].values
            ci = s_data['ci95_half'].values
            npix = s_data['n_distinct_pixels'].values
            
            ax.errorbar(x, y, yerr=ci, label=soil, fmt='-o', color=colors[soil],
                        capsize=5, capthick=1.5, linewidth=2, markersize=6)
            for xi, yi, pxi in zip(x, y, npix):
                if not np.isnan(yi):
                    ax.annotate(f"{pxi}px", (xi, yi), textcoords="offset points",
                                xytext=(0, 6 if soil in ['Loam', 'Clay Loam'] else -10),
                                ha='center', fontsize=6.5, color=colors[soil], alpha=0.85)

        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xticks(range(len(gmia_order)))
        ax.set_xticklabels(gmia_order, fontsize=10)
        ax.set_xlabel("GMIA irrigation equipment class (%)", fontsize=11)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.6, linewidth=1)
        ax.grid(True, linestyle=':', alpha=0.5)
        
    axes[0].set_ylabel("Mean difference $d = \\mathrm{SMAP} - \\mathrm{SSM}_{\\mathrm{OL}}$ [$\\mathrm{m}^3/\\mathrm{m}^3$]", fontsize=11)
    axes[0].legend(title="Dominant soil texture", fontsize=9, title_fontsize=10, loc='lower left')
    plt.suptitle("F3. Mean difference d by GMIA irrigation class, stratified by dominant soil texture and season (95% bootstrap CI & pixel counts)",
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig_path_en = os.path.join(OUT_FIG_DIR_EN, "F3_d_by_gmia_texture_season_en.png")
    plt.savefig(fig_path_en, dpi=300)
    plt.close()
    print(f"  Saved Figure F3 (FR & EN): {fig_path_fr}, {fig_path_en}")


def generate_f4(lat_1d, lon_1d, landmask, df_all):
    """F4: Share of observations rejected by wilting threshold: spatial map and monthly cycle."""
    print("--- Generating Figure F4 & Source CSV ---")
    n_lat, n_lon = len(lat_1d), len(lon_1d)
    
    total_cnt = df_all.groupby('pixel_id').size()
    wilt_cnt = df_all[df_all['qc_lis_wilt_limit']].groupby('pixel_id').size()
    
    rej_rate_grid = np.full((n_lat, n_lon), np.nan, dtype=np.float32)
    csv_map_rows = []
    
    for r in range(n_lat):
        for c in range(n_lon):
            if landmask[r, c] == 1:
                pid = r * n_lon + c
                tot = total_cnt.get(pid, 0)
                wilt = wilt_cnt.get(pid, 0)
                rate = (wilt / tot * 100.0) if tot > 0 else 0.0
                rej_rate_grid[r, c] = rate
                csv_map_rows.append({
                    'pixel_id': pid,
                    'lat': lat_1d[r],
                    'lon': lon_1d[c],
                    'total_valid_smap': tot,
                    'wilt_rejected_obs': wilt,
                    'wilt_rejection_rate_pct': np.round(rate, 2)
                })
                
    df_csv_map = pd.DataFrame(csv_map_rows)
    csv_map_path = os.path.join(OUT_FIG_DIR, "F4_wilt_rejection_map.csv")
    df_csv_map.to_csv(csv_map_path, index=False)
    print(f"  Saved F4 spatial source CSV: {csv_map_path}")
    
    monthly_stats = df_all.groupby(['month', 'pass'])['qc_lis_wilt_limit'].agg(['sum', 'count']).reset_index()
    monthly_stats['rejection_rate_pct'] = (monthly_stats['sum'] / monthly_stats['count']) * 100.0
    
    monthly_global = df_all.groupby('month')['qc_lis_wilt_limit'].agg(['sum', 'count']).reset_index()
    monthly_global['pass'] = 'Global'
    monthly_global['rejection_rate_pct'] = (monthly_global['sum'] / monthly_global['count']) * 100.0
    
    monthly_comb = pd.concat([monthly_stats, monthly_global], ignore_index=True)
    csv_month_path = os.path.join(OUT_FIG_DIR, "F4_wilt_rejection_monthly.csv")
    monthly_comb.to_csv(csv_month_path, index=False)
    print(f"  Saved F4 monthly source CSV: {csv_month_path}")
    
    extent = [lon_1d.min(), lon_1d.max(), lat_1d.min(), lat_1d.max()]
    mask_sea = (landmask == 0)
    rej_m = np.ma.masked_where(mask_sea, rej_rate_grid)
    months = np.arange(1, 13)
    
    am_rates = monthly_stats[monthly_stats['pass'] == 'AM'].sort_values('month')['rejection_rate_pct'].values
    pm_rates = monthly_stats[monthly_stats['pass'] == 'PM'].sort_values('month')['rejection_rate_pct'].values
    glob_rates = monthly_global.sort_values('month')['rejection_rate_pct'].values

    # 1. French Plot
    fig, (ax_map, ax_cycle) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1.2, 1.0]})
    im = ax_map.imshow(rej_m, origin='lower', extent=extent, cmap='YlOrRd', vmin=0, vmax=100)
    ax_map.set_title("A. Taux spatial de rejet par le seuil de flétrissement (%)", fontsize=12, fontweight='bold', pad=10)
    ax_map.set_xlabel("Longitude (°E)", fontsize=10)
    ax_map.set_ylabel("Latitude (°N)", fontsize=10)
    ax_map.grid(True, linestyle='--', alpha=0.4)
    cbar = fig.colorbar(im, ax=ax_map, fraction=0.046, pad=0.04)
    cbar.set_label("Part des observations rejetées (SMAP - θ_wlt < 0.02) [%]", fontsize=9)
    
    month_names_fr = ['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Juin', 'Juil', 'Aoû', 'Sep', 'Oct', 'Nov', 'Déc']
    ax_cycle.plot(months, am_rates, marker='o', linewidth=2.5, color='#1f77b4', label='AM (06h)')
    ax_cycle.plot(months, pm_rates, marker='s', linewidth=2.5, color='#9467bd', label='PM (18h)')
    ax_cycle.plot(months, glob_rates, marker='^', linewidth=2.0, linestyle='--', color='#d62728', label='Moyenne Globale')
    ax_cycle.set_title("B. Cycle mensuel du taux de rejet par flétrissement", fontsize=12, fontweight='bold', pad=10)
    ax_cycle.set_xticks(months)
    ax_cycle.set_xticklabels(month_names_fr, fontsize=10)
    ax_cycle.set_xlabel("Mois de l'année", fontsize=11)
    ax_cycle.set_ylabel("Taux de rejet (%)", fontsize=11)
    ax_cycle.set_ylim(0, 100)
    ax_cycle.grid(True, linestyle=':', alpha=0.5)
    ax_cycle.legend(fontsize=10, loc='upper left')
    
    plt.suptitle("F4. Part des observations rejetées par la contrainte physique de flétrissement (2016-2020)",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig_path_fr = os.path.join(OUT_FIG_DIR, "F4_wilt_rejection.png")
    plt.savefig(fig_path_fr, dpi=300)
    plt.close()
    
    # 2. English Plot (purely descriptive, zero interpretation phrases)
    fig, (ax_map, ax_cycle) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1.2, 1.0]})
    im = ax_map.imshow(rej_m, origin='lower', extent=extent, cmap='YlOrRd', vmin=0, vmax=100)
    ax_map.set_title("A. Spatial rejection rate by wilting threshold (%)", fontsize=12, fontweight='bold', pad=10)
    ax_map.set_xlabel("Longitude (°E)", fontsize=10)
    ax_map.set_ylabel("Latitude (°N)", fontsize=10)
    ax_map.grid(True, linestyle='--', alpha=0.4)
    cbar = fig.colorbar(im, ax=ax_map, fraction=0.046, pad=0.04)
    cbar.set_label("Observations rejected (SMAP - θ_wlt < 0.02) [%]", fontsize=9)
    
    month_names_en = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    ax_cycle.plot(months, am_rates, marker='o', linewidth=2.5, color='#1f77b4', label='AM (06:00)')
    ax_cycle.plot(months, pm_rates, marker='s', linewidth=2.5, color='#9467bd', label='PM (18:00)')
    ax_cycle.plot(months, glob_rates, marker='^', linewidth=2.0, linestyle='--', color='#d62728', label='Global Mean')
    ax_cycle.set_title("B. Monthly cycle of wilting threshold rejection rate", fontsize=12, fontweight='bold', pad=10)
    ax_cycle.set_xticks(months)
    ax_cycle.set_xticklabels(month_names_en, fontsize=10)
    ax_cycle.set_xlabel("Month of year", fontsize=11)
    ax_cycle.set_ylabel("Rejection rate (%)", fontsize=11)
    ax_cycle.set_ylim(0, 100)
    ax_cycle.grid(True, linestyle=':', alpha=0.5)
    ax_cycle.legend(fontsize=10, loc='upper left')
    
    plt.suptitle("F4. Proportion of observations rejected by wilting point proximity threshold (2016-2020)",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig_path_en = os.path.join(OUT_FIG_DIR_EN, "F4_wilt_rejection_en.png")
    plt.savefig(fig_path_en, dpi=300)
    plt.close()
    print(f"  Saved Figure F4 (FR & EN): {fig_path_fr}, {fig_path_en}")


def main():
    os.makedirs(OUT_FIG_DIR, exist_ok=True)
    os.makedirs(OUT_FIG_DIR_EN, exist_ok=True)
    lat_1d, lon_1d, landmask, dom_soil, df_all = load_all_data()
    generate_f1(lat_1d, lon_1d, landmask, df_all)
    generate_f2(lat_1d, lon_1d, landmask, df_all)
    generate_f3(dom_soil, df_all)
    generate_f4(lat_1d, lon_1d, landmask, df_all)
    print("\nDescriptive figures F1-F4 (FR & EN) and CSVs generated successfully.")


if __name__ == "__main__":
    main()
