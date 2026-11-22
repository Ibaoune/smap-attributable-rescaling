# Rapport Final de Synthèse : Phases 1a à 1e (SAR)

## 1. Vue d'Ensemble du Pipeline Réalisé
Le présent rapport documente la finalisation complète des étapes 1a à 1e du projet **smap-attributable-rescaling** (SAR), visant à établir une correction de biais SMAP explicable et attribuable alternative au CDF matching.

- **Domaine spatial** : NorthMor (121 × 181 mailles à 0.05°, 16 030 pixels terrestres, coordonnées vérifiées [30.00°N, 36.00°N] × [-10.00°E, -1.00°E]).
- **Période** : 2016 à 2020 (5 années complètes, 1 827 jours).
- **Total des observations modélisables** : **20,091,609 lignes** réparties sur les 5 années et 2 passages.

## 2. Nombre de Lignes par Année et par Passage

| Année | Passage AM | Passage PM | Total Annuel |
| :---: | :---: | :---: | :---: |
| 2016 | 2,030,375 | 2,012,008 | 4,042,383 |
| 2017 | 2,094,469 | 2,094,030 | 4,188,499 |
| 2018 | 1,960,409 | 1,949,131 | 3,909,540 |
| 2019 | 1,909,153 | 1,906,201 | 3,815,354 |
| 2020 | 2,101,250 | 2,034,583 | 4,135,833 |
| **Total** | **10,095,656** | **9,995,953** | **20,091,609** |

## 3. Statistiques de la Cible $d = \text{SMAP} - \text{SSM}_{\text{OL}}$ par Saison

| Saison | Échantillon ($N$) | Moyenne ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Écart-type ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Médiane ($d$) |
| :--- | :---: | :---: | :---: | :---: |
| **DJF** | 4,951,497 | -0.0583 | 0.0552 | -0.0643 |
| **MAM** | 5,180,356 | -0.0567 | 0.0535 | -0.0583 |
| **JJA** | 4,931,249 | -0.0035 | 0.0440 | +0.0023 |
| **SON** | 5,028,507 | -0.0348 | 0.0522 | -0.0324 |
| **Global** | **20,091,609** | **-0.0386** | **0.0560** | **-0.0357** |

## 4. Statistiques de la Cible par Classe d'Irrigation GMIA

| Classe GMIA | Part échantillon | Moyenne ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Écart-type ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] |
| :--- | :---: | :---: |
| **0%** | 13,482,639 (67.1%) | -0.0362 | 0.0547 |
| **0-10%** | 5,213,447 (25.9%) | -0.0399 | 0.0585 |
| **10-30%** | 733,921 (3.7%) | -0.0543 | 0.0569 |
| **>30%** | 661,602 (3.3%) | -0.0588 | 0.0526 |

## 5. Taux de Rejet des Masques Contrôle Qualité (LIS QC)

| Règle de Contrôle Qualité | Nombre de rejets | Taux de rejet (%) |
| :--- | :---: | :---: |
| `qc_lis_wilt_limit` (Proximité flétrissement (SMAP - θ_wlt < 0.02)) | 9,435,386 | 46.96% |
| `qc_lis_sat_limit` (Proximité saturation (θ_sat - SMAP < 0.02)) | 20,605 | 0.10% |
| `qc_lis_forest` (Couvert forestier (Classes IGBP 1-4)) | 35,393 | 0.18% |
| `qc_lis_rain` (Pluie active (> 0.259 mm/j)) | 3,159,820 | 15.73% |
| `qc_lis_frozen` (Sol gelé (T_soil <= 273.15 K)) | 32,081 | 0.16% |
| `qc_lis_snow` (Neige (SWE > 1 mm ou SCA > 1e-4)) | 193,828 | 0.96% |
| `qc_lis_high_veg` (Forte végétation (shdfac > 0.90)) | 3,699 | 0.02% |
| `qc_lis_vwc_gt_5` (Végétation très dense (VWC > 5.0 kg/m²)) | 0 | 0.00% |
| `qc_lis_combined` (Masque combiné LIS (rejet DA)) | 11,909,408 | 59.28% |

## 6. Validation du Masque Reconstruit vs DAOBS (Enregistrement 2)

Dans le système LIS, les fichiers `LISDAOBS_YYYYMMDDHHMM.1gs4r` comportent deux enregistrements distincts :
1. **Enregistrement 1 (`smobs_unsc`)** : écrit **AVANT** le contrôle qualité `noahmp401_qc_soilmobs` (observation brute reprojectée, filtrée uniquement par le flag SMAP bit 0 et VWC <= 5.0).
2. **Enregistrement 2 (`smobs`)** : écrit **APRÈS** `noahmp401_qc_soilmobs` (pixels rejetés passés à `-9999.0` avant assimilation par EnKF).

- **Nombre de fichiers DAOBS testés** : 122
- **Observations SMAP valides comparées** : 458,464
- **Rejets effectifs observés dans DAOBS (Rec 1 valide, Rec 2 rejeté)** : 221,924 (48.41%)
- **Taux de concordance globale du masque combiné** : **95.48%**
- **Sensibilité de détection des rejets LIS** : **90.67%**

## 7. Contrôle de Qualité et Absence de NaN par Colonne

| Colonne | Type | Statut NaN | Valeurs Manquantes |
| :--- | :---: | :---: | :---: |
| `date` | str | CONFORME (0 NaN) | 0 |
| `year` | int16 | CONFORME (0 NaN) | 0 |
| `month` | int8 | CONFORME (0 NaN) | 0 |
| `season` | str | CONFORME (0 NaN) | 0 |
| `lat` | float32 | CONFORME (0 NaN) | 0 |
| `lon` | float32 | CONFORME (0 NaN) | 0 |
| `pixel_id` | int32 | CONFORME (0 NaN) | 0 |
| `pass` | str | CONFORME (0 NaN) | 0 |
| `pass_is_pm` | int8 | CONFORME (0 NaN) | 0 |
| `spatial_block_id` | str | CONFORME (0 NaN) | 0 |
| `spatial_fold` | int8 | CONFORME (0 NaN) | 0 |
| `split` | str | CONFORME (0 NaN) | 0 |
| `target_d` | float32 | CONFORME (0 NaN) | 0 |
| `smap_sm` | float32 | CONFORME (0 NaN) | 0 |
| `ssm_ol_daily` | float32 | CONFORME (0 NaN) | 0 |
| `elevation` | float32 | CONFORME (0 NaN) | 0 |
| `slope_deg` | float32 | CONFORME (0 NaN) | 0 |
| `slopetype` | int8 | CONFORME (0 NaN) | 0 |
| `vegetation_opacity` | float32 | CONFORME (0 NaN) | 0 |
| `vegetation_water_content` | float32 | CONFORME (0 NaN) | 0 |
| `surface_temperature` | float32 | CONFORME (0 NaN) | 0 |
| `roughness_coefficient` | float32 | ATTENTION (2,535,862 NaN) | 2535862 |
| `static_water_body_fraction` | float32 | CONFORME (0 NaN) | 0 |
| `surface_flag` | uint16 | CONFORME (0 NaN) | 0 |
| `irrigfrac_gmia` | float32 | CONFORME (0 NaN) | 0 |
| `gmia_class` | str | CONFORME (0 NaN) | 0 |
| `precip_mswep` | float32 | CONFORME (0 NaN) | 0 |
| `delta_precip_imerg_mswep` | float32 | CONFORME (0 NaN) | 0 |
| `precip_1d_sum` | float32 | CONFORME (0 NaN) | 0 |
| `precip_3d_sum` | float32 | CONFORME (0 NaN) | 0 |
| `precip_7d_sum` | float32 | CONFORME (0 NaN) | 0 |
| `days_since_rain_gt_1mm` | float32 | CONFORME (0 NaN) | 0 |
| `doy_sin` | float32 | CONFORME (0 NaN) | 0 |
| `doy_cos` | float32 | CONFORME (0 NaN) | 0 |
| `ssm_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `sm_layer2_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `sm_layer3_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `sm_layer4_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `rzsm_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `et_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `surface_runoff_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `baseflow_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `lai_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `ssm_rel_wilt_ref_lag1` | float32 | CONFORME (0 NaN) | 0 |
| `qc_lis_rain` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_frozen` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_forest` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_snow` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_high_veg` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_sat_limit` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_wilt_limit` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_vwc_gt_5` | bool | CONFORME (0 NaN) | 0 |
| `qc_lis_combined` | bool | CONFORME (0 NaN) | 0 |
| `texture_frac_1` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_2` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_3` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_4` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_5` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_6` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_7` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_8` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_9` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_10` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_11` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_12` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_13` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_14` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_15` | float32 | CONFORME (0 NaN) | 0 |
| `texture_frac_16` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_1` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_2` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_3` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_4` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_5` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_6` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_7` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_8` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_9` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_10` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_11` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_12` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_13` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_14` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_15` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_16` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_17` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_18` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_19` | float32 | CONFORME (0 NaN) | 0 |
| `landcover_frac_20` | float32 | CONFORME (0 NaN) | 0 |

## 8. Points 'À Vérifier' et Choix Méthodologiques
1. **Heure du passage vs Moyenne journalière de l'OL (Point ouvert, Décision 1)** : La cible $d$ utilise actuellement la moyenne journalière `SoilMoist_tavg` de l'OL. L'impact par rapport à un état OL extrait à l'heure précise du passage (06h / 18h) sera quantifié dans l'analyse de sensibilité.
2. **Pluie antécédente IMERG demi-horaire (Décision 7)** : Les cumuls antécédents `precip_1d_sum`, `precip_3d_sum`, `precip_7d_sum` et `days_since_rain_gt_1mm` sont strictement calculés sur les jours $[t-n, t-1]$, éliminant tout risque de fuite de la pluie contemporaine du jour $t$.
3. **Conservation MSWEP 0.10° vers NorthMor 0.05°** : L'affectation directe conservative au plus proche voisin préserve strictement les volumes totaux sur les bassins hydrographiques.
4. **Conservation des composantes cibles comme non-prédicteurs** : `smap_sm` et `ssm_ol_daily` figurent dans la table pour traçabilité et décomposition physique, mais portent le statut strict de non-prédicteurs.
