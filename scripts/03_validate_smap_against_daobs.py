#!/usr/bin/env python3
"""
Step 1b Validation: Compare SAR raw reprojected SMAP against unscaled DAOBS binary files.
Sampling: Exactly 2 observation files per month and per pass (AM and PM) over 2016-2020 (60 months x 2 passes x 2 = up to 240 files).
Produces dynamic report with conclusion generated strictly from the numbers.
"""

import os
import sys
import glob
import re
import datetime
from collections import defaultdict
import numpy as np
import netCDF4

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from sar.regrid import latlon_to_ease2_m09, regrid_smap_to_northmor
from sar.smap import open_smap_h5, read_smap_pass_data

DAOBS_ROOT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/DA_nocdf_noirr_2016_2020/output"
SMAP_RAW_DIR = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/observations/SMAP_SPL3SMP_E/raw"
LIS_INPUT = "/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc"
REPORT_FILE = os.path.join(os.path.dirname(__file__), '..', 'reports', 'validation_1b_daobs_report.md')


def read_daobs_file(filepath: str, n_lat: int = 121, n_lon: int = 181) -> np.ndarray:
    """
    Read first record (smobs_unsc) of LIS sequential binary DAOBS file.
    Shape in Fortran is (181, 121) column-major, returned as (121, 181).
    """
    try:
        with open(filepath, 'rb') as f:
            rec_hdr = np.fromfile(f, dtype='>u4', count=1)
            if len(rec_hdr) == 0:
                return None
            data = np.fromfile(f, dtype='>f4', count=n_lat * n_lon)
            if len(data) != n_lat * n_lon:
                return None
            grid = data.reshape((n_lat, n_lon))
            return grid
    except Exception as e:
        return None


def get_grid_mapping():
    ds = netCDF4.Dataset(LIS_INPUT, 'r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    ds.close()
    return latlon_to_ease2_m09(lat, lon)


def find_smap_raw_file(date_str: str) -> str:
    matches = glob.glob(os.path.join(SMAP_RAW_DIR, f"SMAP_L3_SM_P_E_{date_str}_*.h5"))
    return matches[0] if matches else None


def main():
    print("Loading grid mapping...")
    row_idx, col_idx = get_grid_mapping()

    print("Locating DAOBS files...")
    daobs_pattern = os.path.join(DAOBS_ROOT, "*", "DAOBS", "*", "LISDAOBS_*.1gs4r")
    daobs_files = sorted(glob.glob(daobs_pattern))
    print(f"Total DAOBS files located: {len(daobs_files)}")

    # Index files by (year, month, pass)
    files_by_group = defaultdict(list)
    for fpath in daobs_files:
        fname = os.path.basename(fpath)
        match = re.search(r'LISDAOBS_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})', fname)
        if not match:
            continue
        yyyy, mm, dd, hh, mn = match.groups()
        hour = int(hh)
        pass_name = 'AM' if hour < 12 else 'PM'
        group_key = (int(yyyy), int(mm), pass_name)
        files_by_group[group_key].append(fpath)

    # Sampling: 2 files per month and per pass over 2016-2020
    years = [2016, 2017, 2018, 2019, 2020]
    target_groups = []
    for y in years:
        for m in range(1, 13):
            for p in ['AM', 'PM']:
                target_groups.append((y, m, p))

    print(f"Total target monthly groups to sample: {len(target_groups)} (2 files each -> up to {len(target_groups)*2} files)")

    total_pixels_compared = 0
    exact_matches_count = 0
    diff_1e6_count = 0
    all_diffs = []
    validated_files = []
    smap_cache = {}

    for (y, m, p) in target_groups:
        candidate_files = files_by_group.get((y, m, p), [])
        valid_in_group = 0
        
        for fpath in candidate_files:
            fname = os.path.basename(fpath)
            match = re.search(r'LISDAOBS_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})', fname)
            yyyy, mm, dd, hh, mn = match.groups()
            date_str = f"{yyyy}{mm}{dd}"

            daobs_data = read_daobs_file(fpath)
            if daobs_data is None:
                continue

            valid_mask = daobs_data != -9999.0
            n_valid = int(np.sum(valid_mask))
            if n_valid == 0:
                continue

            cache_key = (date_str, p)
            if cache_key not in smap_cache:
                smap_path = find_smap_raw_file(date_str)
                if smap_path is None or not os.path.exists(smap_path):
                    continue
                try:
                    ds_smap = open_smap_h5(smap_path)
                    smap_pass_data = read_smap_pass_data(ds_smap, p)
                    ds_smap.close()

                    sm_raw = smap_pass_data['soil_moisture']
                    qa_raw = smap_pass_data['retrieval_qual_flag']
                    vwc_raw = smap_pass_data['vegetation_water_content']

                    # Replicate exact LIS QC from read_NASASMAPsm.F90:
                    # 1. vwc <= 5.0
                    # 2. retrieval_qual_flag bit 0 == 0
                    sm_filtered = sm_raw.copy()
                    sm_filtered[vwc_raw > 5.0] = np.nan
                    sm_filtered[(qa_raw.astype(np.uint16) & 1) != 0] = np.nan

                    pred_regrid = regrid_smap_to_northmor(sm_filtered, row_idx, col_idx)
                    smap_cache[cache_key] = pred_regrid
                except Exception as e:
                    print(f"Error reading SMAP {date_str} {p}: {e}")
                    continue

            pred_regrid = smap_cache[cache_key]

            # Compare on all valid DAOBS pixels
            daobs_vals = daobs_data[valid_mask]
            pred_vals = pred_regrid[valid_mask]

            # Differences
            # If pred_vals has NaN where daobs has valid, difference is NaN/mismatch
            diffs = np.abs(pred_vals - daobs_vals)
            # If any NaN in pred_vals:
            nan_in_pred = np.isnan(pred_vals)
            if np.any(nan_in_pred):
                print(f"Warning: {np.sum(nan_in_pred)} NaNs in prediction where DAOBS was valid in {fname}")

            n_exact = int(np.sum(diffs == 0.0))
            n_1e6 = int(np.sum(diffs < 1e-6))
            max_d = float(np.nanmax(diffs)) if len(diffs) > 0 else 0.0

            total_pixels_compared += n_valid
            exact_matches_count += n_exact
            diff_1e6_count += n_1e6
            all_diffs.extend(diffs[~np.isnan(diffs)])

            valid_in_group += 1
            validated_files.append({
                'file': fname,
                'date': f"{yyyy}-{mm}-{dd}",
                'month': f"{yyyy}-{mm}",
                'pass': p,
                'valid_pixels': n_valid,
                'exact_matches': n_exact,
                'max_diff': max_d
            })

            if valid_in_group >= 2:
                break

    all_diffs = np.array(all_diffs)
    exact_rate = 100.0 * exact_matches_count / total_pixels_compared if total_pixels_compared > 0 else 0.0
    rate_1e6 = 100.0 * diff_1e6_count / total_pixels_compared if total_pixels_compared > 0 else 0.0
    max_diff = float(np.max(all_diffs)) if len(all_diffs) > 0 else 0.0
    mean_diff = float(np.mean(all_diffs)) if len(all_diffs) > 0 else 0.0

    print("\n================ VALIDATION RESULTS 1b ================")
    print(f"Total DAOBS files evaluated: {len(validated_files)}")
    print(f"Total valid observation pixels compared: {total_pixels_compared:,}")
    print(f"Exact matches (diff == 0.0): {exact_matches_count:,} ({exact_rate:.4f}%)")
    print(f"Matches within float tolerance (diff < 1e-6): {diff_1e6_count:,} ({rate_1e6:.4f}%)")
    print(f"Maximum difference: {max_diff:.8e}")
    print(f"Mean absolute difference: {mean_diff:.8e}")
    print("=======================================================\n")

    # Dynamic conclusion generation
    if total_pixels_compared > 0 and exact_rate == 100.0:
        conclusion_text = (
            f"La validation exhaustive sur **{len(validated_files)} fichiers** (échantillon systématique de 2 fichiers par mois et par passage AM/PM sur l'ensemble de la période 2016–2020) "
            f"confirme une concordance numérique absolue. Sur un total de **{total_pixels_compared:,} pixels-observations comparés**, "
            f"le nombre d'égalités strictes (`diff == 0.0`) est exactement de **{exact_matches_count:,} ({exact_rate:.4f}%)**, "
            f"avec un écart moyen de `{mean_diff:.2e}` et un écart maximal de `{max_diff:.2e}` $\\mathrm{{m}}^3/\\mathrm{{m}}^3$. "
            f"La chaîne de reprojection au plus proche voisin de SAR reproduit donc sans aucune approximation la routine Fortran `neighbor_interp_input_withgrid.F90` de LIS."
        )
    else:
        conclusion_text = (
            f"Sur un total de {total_pixels_compared:,} pixels-observations comparés à travers {len(validated_files)} fichiers, "
            f"le taux d'accord à 1e-6 est de {rate_1e6:.4f}% ({diff_1e6_count:,} pixels), "
            f"avec un écart maximal de {max_diff:.4e} et un écart moyen de {mean_diff:.4e}."
        )

    # Generate Markdown Report
    os.makedirs(os.path.dirname(REPORT_FILE), exist_ok=True)
    with open(REPORT_FILE, 'w') as rf:
        rf.write("# Rapport de Validation 1b : Reprojection SMAP vs DAOBS LIS\n\n")
        rf.write("## 1. Contexte et Méthodologie\n")
        rf.write("Pour valider rigoureusement la chaîne de reprojection au plus proche voisin (`src/sar/regrid.py`), ")
        rf.write("les observations SMAP Enhanced brutes (`SPL3SMP_E`) reprojectées ont été comparées ")
        rf.write("cellule par cellule aux observations non mises à l'échelle (`smobs_unsc`, Enregistrement 1) lues directement ")
        rf.write("dans les fichiers binaires séquentiels `DAOBS` issus de la simulation LIS `DA_nocdf_noirr_2016_2020`.\n\n")
        rf.write("Le protocole d'échantillonnage sélectionne systématiquement **2 fichiers par mois et par passage (AM et PM)** ")
        rf.write("sur l'intégralité des 5 années (2016 à 2020, soit 60 mois x 2 passages x 2 = jusqu'à 240 fichiers).\n\n")

        rf.write("## 2. Synthèse Globale des Écarts\n")
        rf.write(f"- **Fichiers DAOBS échantillonnés avec observations** : {len(validated_files)}\n")
        rf.write(f"- **Nombre total de pixels-jours comparés** : {total_pixels_compared:,}\n")
        rf.write(f"- **Taux d'égalité stricte (`diff == 0.0`)** : **{exact_rate:.4f}%** ({exact_matches_count:,} / {total_pixels_compared:,})\n")
        rf.write(f"- **Taux d'égalité à précision machine (`diff < 1e-6`)** : **{rate_1e6:.4f}%** ({diff_1e6_count:,} / {total_pixels_compared:,})\n")
        rf.write(f"- **Écart absolu maximal** : `{max_diff:.8e}` $\\mathrm{{m}}^3/\\mathrm{{m}}^3$\n")
        rf.write(f"- **Écart absolu moyen** : `{mean_diff:.8e}` $\\mathrm{{m}}^3/\\mathrm{{m}}^3$\n\n")

        rf.write("## 3. Répartition des Écarts\n")
        rf.write("| Seuil d'écart | Nombre de pixels | Pourcentage cumulé |\n")
        rf.write("| :--- | :---: | :---: |\n")
        rf.write(f"| `diff == 0.0` (stricte égalité) | {exact_matches_count:,} | {exact_rate:.4f}% |\n")
        rf.write(f"| `diff < 1e-7` | {np.sum(all_diffs < 1e-7):,} | {100.0 * np.sum(all_diffs < 1e-7) / total_pixels_compared:.4f}% |\n")
        rf.write(f"| `diff < 1e-6` | {diff_1e6_count:,} | {rate_1e6:.4f}% |\n")
        rf.write(f"| `diff >= 1e-6` | {np.sum(all_diffs >= 1e-6):,} | {100.0 * np.sum(all_diffs >= 1e-6) / total_pixels_compared:.4f}% |\n\n")

        rf.write("## 4. Échantillon par Année et Mois (Extrait)\n")
        rf.write("| Fichier DAOBS | Date | Mois | Passage | Pixels observés | Égalités exactes | Écart max |\n")
        rf.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for item in validated_files[:40]:
            rf.write(f"| `{item['file']}` | {item['date']} | {item['month']} | {item['pass']} | {item['valid_pixels']} | {item['exact_matches']} | {item['max_diff']:.2e} |\n")
        if len(validated_files) > 40:
            rf.write(f"| ... ({len(validated_files) - 40} autres fichiers analysés) | ... | ... | ... | ... | ... | ... |\n")

        rf.write("\n## 5. Conclusion\n")
        rf.write(f"{conclusion_text}\n")

    print(f"Updated validation report written to: {REPORT_FILE}")

if __name__ == "__main__":
    main()
