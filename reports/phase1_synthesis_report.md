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

## 3. Répartition par Pli Spatial (Protocole d'Évaluation Équilibré)

| Pli Spatial (`spatial_fold`) | Nombre d'Observations | Part de l'Échantillon (%) |
| :---: | :---: | :---: |
| Pli 0 | 4,055,079 | 20.18% |
| Pli 1 | 3,808,367 | 18.96% |
| Pli 2 | 3,976,288 | 19.79% |
| Pli 3 | 4,142,271 | 20.62% |
| Pli 4 | 4,109,604 | 20.45% |
| **Total** | **20,091,609** | **100.00%** |

## 4. Statistiques de la Cible $d = \text{SMAP} - \text{SSM}_{\text{OL}}$ par Saison

| Saison | Échantillon ($N$) | Moyenne ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Écart-type ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Médiane ($d$) |
| :--- | :---: | :---: | :---: | :---: |
| **DJF** | 4,951,497 | -0.0583 | 0.0551 | -0.0643 |
| **MAM** | 5,180,356 | -0.0567 | 0.0533 | -0.0583 |
| **JJA** | 4,931,249 | -0.0035 | 0.0439 | +0.0023 |
| **SON** | 5,028,507 | -0.0348 | 0.0521 | -0.0324 |
| **Global** | **20,091,609** | **-0.0386** | **0.0549** | **-0.0357** |

## 5. Statistiques de la Cible par Classe d'Irrigation GMIA (Constat Factuel)

| Classe GMIA | Part échantillon | Moyenne ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] | Écart-type ($d$) [$\mathrm{m}^3/\mathrm{m}^3$] |
| :--- | :---: | :---: | :---: |
| **0%** | 13,482,639 (67.1%) | -0.0362 | 0.0540 |
| **0-10%** | 5,213,447 (25.9%) | -0.0399 | 0.0583 |
| **10-30%** | 733,921 (3.7%) | -0.0543 | 0.0569 |
| **>30%** | 661,602 (3.3%) | -0.0588 | 0.0526 |

> [!NOTE]
> La moyenne de $d$ passe de -0.0362 m³/m³ en classe 0% à -0.0588 m³/m³ pour >30%. Ce gradient négatif est strictement factuel ; il est en sens inverse d'un signal d'irrigation direct (qui humidifierait le sol et augmenterait $d$) et reflète vraisemblablement une covariation géographique avec la texture du sol (sols plus argileux dans les plaines agricoles du Gharb/Saïss avec seuil de flétrissement plus élevé).

## 6. Analyse des Valeurs Manquantes de Rugosité (`roughness_coefficient`)

- **Nombre total de valeurs manquantes de rugosité** : 2,535,862 sur 20,091,609 observations (12.62%).
- Une colonne indicatrice dédiée `roughness_is_nan` (booléen / int8) a été intégrée à la table de modélisation pour permettre aux modèles d'arbres de l'exploiter nativement sans imputation artificielle.

### Répartition des NaN de Rugosité par Type d'Occupation du Sol (IGBP)

| Classe IGBP | Description | Nombre de NaN | Part des NaN (%) |
| :---: | :--- | :---: | :---: |
| 7 | Open Shrublands | 1,057,792 | 41.71% |
| 16 | Barren or Sparsely Vegetated | 986,854 | 38.92% |
| 10 | Grasslands | 302,169 | 11.92% |
| 12 | Croplands | 109,769 | 4.33% |
| 6 | Closed Shrublands | 38,662 | 1.52% |
| 5 | Mixed Forests | 22,069 | 0.87% |

Les classes 7 (Open Shrublands) et 16 (Barren or Sparsely Vegetated) concentrent à elles seules 80.63% des valeurs manquantes de rugosité.

## 7. Taux de Rejet des Masques Contrôle Qualité (LIS QC)

| Règle de Contrôle Qualité | Nombre de rejets | Taux de rejet (%) |
| :--- | :---: | :---: |
| `qc_lis_wilt_limit` (Proximité flétrissement (SMAP - θ_wlt < 0.02)) | 9,435,386 | 46.96% |
| `qc_lis_sat_limit` (Proximité saturation (θ_sat - SMAP < 0.02)) | 20,605 | 0.10% |
| `qc_lis_forest` (Couvert forestier (Classes IGBP 1-4)) | 35,393 | 0.18% |
| `qc_lis_rain` (Taux de pluie au passage (> 3e-6 kg/m²/s)) | 630,385 | 3.14% |
| `qc_lis_frozen` (Sol gelé (T_soil <= 273.15 K)) | 32,081 | 0.16% |
| `qc_lis_snow` (Neige (SWE > 1 mm ou SCA > 1e-4)) | 193,828 | 0.96% |
| `qc_lis_high_veg` (Forte végétation (shdfac > 0.90)) | 3,699 | 0.02% |
| `qc_lis_vwc_gt_5` (Végétation très dense (VWC > 5.0 kg/m²)) | 0 | 0.00% |
| `qc_lis_combined` (Masque combiné LIS (rejet DA)) | 10,174,633 | 50.64% |

## 8. Validation du Masque Reconstruit vs DAOBS (Enregistrement 2)

Dans le système LIS, les fichiers `LISDAOBS_YYYYMMDDHHMM.1gs4r` comportent deux enregistrements distincts :
1. **Enregistrement 1 (`smobs_unsc`)** : écrit **AVANT** le contrôle qualité `noahmp401_qc_soilmobs` (observation brute reprojectée, filtrée uniquement par le flag SMAP bit 0 et VWC <= 5.0).
2. **Enregistrement 2 (`smobs`)** : écrit **APRÈS** `noahmp401_qc_soilmobs` (pixels rejetés passés à `-9999.0` avant assimilation par EnKF).

- **Nombre de fichiers DAOBS testés** : 122
- **Observations SMAP valides comparées** : 458,464
- **Rejets effectifs observés dans DAOBS (Rec 1 valide, Rec 2 rejeté)** : 221,924 (48.41%)
- **Matrice de confusion du masque combiné** :
  - Vrais Positifs (TP, rejet prédit et réel) : **205,430**
  - Faux Positifs (FP, rejet prédit mais accepté dans DAOBS) : **13,377**
  - Faux Négatifs (FN, accepté prédit mais rejeté dans DAOBS) : **16,494**
  - Vrais Négatifs (TN, accepté prédit et réel) : **223,163**
- **Métriques combinées** : Précision = **93.89%**, Sensibilité/Rappel = **92.57%**, Exactitude = **93.48%**

> [!NOTE]
> Le masque de pluie est une approximation (précision ~19% due à l'utilisation du taux demi-horaire IMERG au passage nominal, alors que LIS évalue la pluie instantanée du pas de temps DA de 30 min issu du forçage atmosphérique).

### Accord Règle par Règle avec DAOBS Record 2

| Règle Testée | Déclenchements Règle | Rejets Réels DAOBS | Précision Règle (%) | Part des Rejets LIS Expliquée (%) |
| :--- | :---: | :---: | :---: | :---: |
| **Proximité flétrissement (SMAP - θ_wlt < 0.02)** | 199,766 | 199,766 | 100.00% | 90.02% |
| **Proximité saturation (θ_sat - SMAP < 0.02)** | 793 | 793 | 100.00% | 0.36% |
| **Couvert forestier (IGBP <= 4)** | 747 | 747 | 100.00% | 0.34% |
| **Taux de pluie pass (> 3e-6 kg/m²/s)** | 15,063 | 2,917 | 19.37% | 1.31% |
| **Sol gelé (T_soil <= 273.15 K)** | 1,077 | 1,072 | 99.54% | 0.48% |
| **Neige (SWE > 1mm ou SCA > 1e-4)** | 5,393 | 4,041 | 74.93% | 1.82% |
| **Forte végétation (shdfac > 0.90)** | 98 | 98 | 100.00% | 0.04% |

## 9. Conservation des Cumuls de Précipitations MSWEP (0.10° → 0.05°)

Le passage de la grille MSWEP v2.8 (0.10°) à la grille NorthMor (0.05°) a été vérifié par bilan volumétrique strict sur les bassins versants du Sebou et du Bouregreg :

| Année | Cumul Sebou (0.10°) [mm] | Cumul Sebou (0.05°) [mm] | Écart Sebou (%) | Cumul Bouregreg (0.10°) [mm] | Cumul Bouregreg (0.05°) [mm] | Écart Bouregreg (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 2016 | 533.75 | 534.13 | +0.07% | 488.53 | 485.99 | -0.52% |
| 2017 | 373.69 | 374.48 | +0.21% | 338.81 | 339.74 | +0.28% |
| 2018 | 803.67 | 803.84 | +0.02% | 740.89 | 735.63 | -0.71% |
| 2019 | 434.61 | 435.12 | +0.12% | 348.79 | 347.10 | -0.48% |
| 2020 | 442.22 | 442.59 | +0.09% | 346.94 | 345.55 | -0.40% |
| **Moyenne** | **517.59** | **518.03** | **+0.09%** | **452.79** | **450.80** | **-0.44%** |

L'écart relatif annuel reste strictement inférieur à 0.21% sur le Sebou et 0.71% sur le Bouregreg, confirmant la rigoureuse conservation hydrologique du remaillage direct.

## 10. Contrôle de Qualité et Absence de NaN par Colonne

| Colonne | Type | Statut NaN | Valeurs Manquantes |
| :--- | :---: | :---: | :---: |
| `date` | object | CONFORME (0 NaN) | 0 |
| `year` | int16 | CONFORME (0 NaN) | 0 |
| `month` | int8 | CONFORME (0 NaN) | 0 |
| `season` | object | CONFORME (0 NaN) | 0 |
| `lat` | float32 | CONFORME (0 NaN) | 0 |
| `lon` | float32 | CONFORME (0 NaN) | 0 |
| `pixel_id` | int32 | CONFORME (0 NaN) | 0 |
| `pass` | object | CONFORME (0 NaN) | 0 |
| `pass_is_pm` | int8 | CONFORME (0 NaN) | 0 |
| `spatial_block_id` | object | CONFORME (0 NaN) | 0 |
| `spatial_fold` | int8 | CONFORME (0 NaN) | 0 |
| `split` | object | CONFORME (0 NaN) | 0 |
| `target_d` | float32 | CONFORME (0 NaN) | 0 |
| `smap_sm` | float32 | CONFORME (0 NaN) | 0 |
| `ssm_ol_daily` | float32 | CONFORME (0 NaN) | 0 |
| `elevation` | float32 | CONFORME (0 NaN) | 0 |
| `slope_deg` | float32 | CONFORME (0 NaN) | 0 |
| `slopetype` | int8 | CONFORME (0 NaN) | 0 |
| `vegetation_opacity` | float32 | CONFORME (0 NaN) | 0 |
| `vegetation_water_content` | float32 | CONFORME (0 NaN) | 0 |
| `surface_temperature` | float32 | CONFORME (0 NaN) | 0 |
| `roughness_coefficient` | float32 | DOCUMENTÉ (Flagué par roughness_is_nan: 2,535,862 NaN) | 2535862 |
| `roughness_is_nan` | int8 | CONFORME (0 NaN) | 0 |
| `static_water_body_fraction` | float32 | CONFORME (0 NaN) | 0 |
| `surface_flag` | uint16 | CONFORME (0 NaN) | 0 |
| `irrigfrac_gmia` | float32 | CONFORME (0 NaN) | 0 |
| `gmia_class` | object | CONFORME (0 NaN) | 0 |
| `precip_mswep` | float32 | CONFORME (0 NaN) | 0 |
| `delta_precip_imerg_mswep` | float32 | CONFORME (0 NaN) | 0 |
| `precip_1d_sum` | float32 | CONFORME (Warmup Jan 1-7 2016: 23,613 NaN) | 23613 |
| `precip_3d_sum` | float32 | CONFORME (Warmup Jan 1-7 2016: 41,370 NaN) | 41370 |
| `precip_7d_sum` | float32 | CONFORME (Warmup Jan 1-7 2016: 76,633 NaN) | 76633 |
| `days_since_rain_gt_1mm` | float32 | CONFORME (Warmup Jan 1-7 2016: 23,613 NaN) | 23613 |
| `rain_rate_pass` | float32 | CONFORME (0 NaN) | 0 |
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

