# Rapport de Phase 2c : Correction Attribuable LightGBM et TreeSHAP Exact

**Date d'exécution :** 2026-11-23 18:45:19  
**Période d'entraînement :** 2016-01-01 au 2018-12-31 (scindée en 5 plis spatiaux de 1°)  
**Période d'évaluation :** 2019-01-01 au 2019-12-31 (année 2019 scellée au préalable)  
**Année 2020 :** Strictement scellée (aucune lecture, aucune inférence)  
**Population d'entraînement :** Échantillon stratifié d'au plus 100 observations par pixel sans remise parmi les observations non-wilt QC-valides sur 2016–2018 (min : 100, médiane : 100.0, 0 pixels sous 100 obs sur 14,250 pixels, 1,425,000 lignes au total).  
**Algorithme principal :** LightGBM (300 arbres, learning_rate=0,05, num_leaves=63)  
**Gestion des valeurs manquantes :** Traitement natif des NaN par LightGBM pour la rugosité (`roughness_coefficient`), sans imputation par zéro, complété par l'indicateur binaire `roughness_is_nan`.  
**Méthode d'attribution :** TreeSHAP EXACT natif (`pred_contrib=True`), aucune approximation de Saabas.  
**Constante de base :** $\phi_0 = -0.0430$ m³/m³ (déviation de l'additivité Lundberg : $\le 10^{-15}$).

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
| Pli 0 | 1,132,300 | 292,700 | 0.0212 m³/m³ | 0.9292 | -0.0020 m³/m³ |
| Pli 1 | 1,142,000 | 283,000 | 0.0214 m³/m³ | 0.9339 | -0.0014 m³/m³ |
| Pli 2 | 1,150,800 | 274,200 | 0.0220 m³/m³ | 0.9114 | -0.0007 m³/m³ |
| Pli 3 | 1,135,400 | 289,600 | 0.0206 m³/m³ | 0.9246 | -0.0002 m³/m³ |
| Pli 4 | 1,139,500 | 285,500 | 0.0224 m³/m³ | 0.9095 | +0.0025 m³/m³ |
| **Moyenne** | — | — | **0.0215 m³/m³** | **0.9217** | — |

---

## 3. Stabilité des attributions SHAP (Groupes, Intra-groupe et Variables)

### 3.1 Stabilité au niveau des 4 GROUPES et corrélation de rang intra-groupe

| Groupe | Contribution Absolue Moyenne | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | Rang Moyen | Spearman intra-groupe moyen |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `context_representation` | 0.0431 m³/m³ | 1 | 1 | 1 | 1 | 1 | 1.00 | 0.8800 |
| `removable` | 0.0396 m³/m³ | 2 | 2 | 2 | 2 | 2 | 2.00 | 0.9709 |
| `context_process` | 0.0148 m³/m³ | 3 | 3 | 3 | 3 | 3 | 3.00 | 0.8476 |
| `retained` | 0.0072 m³/m³ | 4 | 4 | 4 | 4 | 4 | 4.00 | 0.9786 |
| **Moyenne groupes** | — | — | — | — | — | — | — | **1.0000** |

### 3.2 Rangs des 15 variables principales à travers les plis

| Rang | Variable | Rang Moyen | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | SHAP Absolu Moyen |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | `ssm_lag1` | 1.00 | 1 | 1 | 1 | 1 | 1 | 0.0296 |
| 2 | `roughness_coefficient` | 2.00 | 2 | 2 | 2 | 2 | 2 | 0.0126 |
| 3 | `vegetation_water_content` | 3.40 | 3 | 3 | 3 | 4 | 4 | 0.0062 |
| 4 | `ssm_rel_wilt_ref_lag1` | 4.20 | 5 | 5 | 5 | 3 | 3 | 0.0059 |
| 5 | `doy_sin` | 4.40 | 4 | 4 | 4 | 5 | 5 | 0.0056 |
| 6 | `pass_is_pm` | 6.00 | 6 | 6 | 6 | 6 | 6 | 0.0045 |
| 7 | `vegetation_opacity` | 7.80 | 7 | 7 | 8 | 7 | 10 | 0.0036 |
| 8 | `surface_temperature` | 8.60 | 9 | 8 | 10 | 8 | 8 | 0.0034 |
| 9 | `lai_lag1` | 9.40 | 10 | 9 | 9 | 12 | 7 | 0.0032 |
| 10 | `sm_layer4_lag1` | 9.80 | 11 | 11 | 7 | 9 | 11 | 0.0030 |
| 11 | `sm_layer3_lag1` | 9.80 | 8 | 10 | 12 | 10 | 9 | 0.0031 |
| 12 | `et_lag1` | 12.00 | 12 | 12 | 11 | 11 | 14 | 0.0024 |
| 13 | `precip_mswep` | 13.20 | 13 | 13 | 14 | 13 | 13 | 0.0021 |
| 14 | `doy_cos` | 14.40 | 15 | 14 | 13 | 15 | 15 | 0.0019 |
| 15 | `precip_3d_sum` | 16.00 | 16 | 15 | 17 | 16 | 16 | 0.0016 |

### 3.3 Matrice de corrélation de Spearman entre plis

| Paire de plis | Spearman $\rho$ (4 Groupes) | Spearman $\rho$ (Top 15) | Spearman $\rho$ (Toutes les 66 variables) |
| :--- | :--- | :--- | :--- |
| Pli 0 vs Pli 1 | 1.0000 | 0.9893 (p=3.00e-12) | 0.9846 (p=3.65e-50) |
| Pli 0 vs Pli 2 | 1.0000 | 0.9321 (p=4.24e-07) | 0.9863 (p=7.42e-52) |
| Pli 0 vs Pli 3 | 1.0000 | 0.9643 (p=7.07e-09) | 0.9893 (p=3.63e-55) |
| Pli 0 vs Pli 4 | 1.0000 | 0.9500 (p=6.09e-08) | 0.9792 (p=4.76e-46) |
| Pli 1 vs Pli 2 | 1.0000 | 0.9500 (p=6.09e-08) | 0.9798 (p=1.91e-46) |
| Pli 1 vs Pli 3 | 1.0000 | 0.9643 (p=7.07e-09) | 0.9909 (p=1.93e-57) |
| Pli 1 vs Pli 4 | 1.0000 | 0.9607 (p=1.30e-08) | 0.9728 (p=2.37e-42) |
| Pli 2 vs Pli 3 | 1.0000 | 0.9464 (p=9.45e-08) | 0.9828 (p=1.08e-48) |
| Pli 2 vs Pli 4 | 1.0000 | 0.9071 (p=3.06e-06) | 0.9727 (p=2.48e-42) |
| Pli 3 vs Pli 4 | 1.0000 | 0.9214 (p=1.07e-06) | 0.9757 (p=6.81e-44) |
| **Moyenne** | **1.0000** | **0.9486** | **0.9814** |

---

## 4. Variantes de correction et décomposition saisonnière (2019)

Les 3 variantes de correction appliquent la décomposition additive :
- **V1 (removable) :** $\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$
- **V2 (removable + context_rep) :** $\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$
- **V3 (removable + all context) :** $\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$

### Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ et biais résiduel $\overline{d_{\mathrm{res}}}$ par saison

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{\mathrm{res}}$ | V2 : Correction | V2 : Biais $d_{\mathrm{res}}$ | V3 : Correction | V3 : Biais $d_{\mathrm{res}}$ | B2-train : Correction | B2-train : Biais $d_{\mathrm{res}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| DJF | 994,881 | +0.0406 m³/m³ | -0.0130 m³/m³ | +0.0569 m³/m³ | +0.0033 m³/m³ | +0.0561 m³/m³ | +0.0024 m³/m³ | +0.0709 m³/m³ | +0.0162 m³/m³ |
| MAM | 1,029,610 | +0.0449 m³/m³ | -0.0139 m³/m³ | +0.0606 m³/m³ | +0.0018 m³/m³ | +0.0622 m³/m³ | +0.0034 m³/m³ | +0.0652 m³/m³ | +0.0046 m³/m³ |
| JJA | 638,878 | +0.0451 m³/m³ | +0.0399 m³/m³ | +0.0053 m³/m³ | +0.0001 m³/m³ | +0.0027 m³/m³ | -0.0025 m³/m³ | +0.0061 m³/m³ | +0.0007 m³/m³ |
| SON | 997,430 | +0.0470 m³/m³ | +0.0207 m³/m³ | +0.0352 m³/m³ | +0.0089 m³/m³ | +0.0310 m³/m³ | +0.0047 m³/m³ | +0.0441 m³/m³ | +0.0170 m³/m³ |
| Annual | 3,660,799 | +0.0443 m³/m³ | +0.0052 m³/m³ | +0.0430 m³/m³ | +0.0039 m³/m³ | +0.0416 m³/m³ | +0.0025 m³/m³ | +0.0505 m³/m³ | +0.0103 m³/m³ |

---

## 5. Protocole 2a complet : Comparaison avec B0, B1, B2 et B2-train

Le flétrissement ($\mathrm{SMAP}^* \ge \theta_{\mathrm{wlt}} + 0{,}02$) est recalculé sur l'observation corrigée pour chaque référence.

### 5.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,912,072 (52.2%) | 1,748,727 (47.8%) | -0.0265 m³/m³ | 0.0597 m³/m³ | 0.0535 m³/m³ | 0.6278 | 1.0000 | +0.0109 [-0.0104, +0.0187] (interprétable) |
| B1_Climatology | 3,133,497 (85.6%) | 527,302 (14.4%) | +0.0132 m³/m³ | 0.0427 m³/m³ | 0.0406 m³/m³ | 0.7605 | 1.0000 | +0.0046 [+0.0003, +0.0093] (interprétable) |
| B2_CDF_Legacy | 2,622,380 (71.6%) | 652,539 (17.8%) | +0.0086 m³/m³ | 0.0354 m³/m³ | 0.0343 m³/m³ | 0.8040 | 0.8009 | +0.0010 [-0.0004, +0.0031] (interprétable) |
| B2_train_CDF | 2,787,534 (76.1%) | 596,740 (16.3%) | +0.0168 m³/m³ | 0.0433 m³/m³ | 0.0399 m³/m³ | 0.7451 | 0.7801 | +0.0026 [+0.0003, +0.0048] (interprétable) |
| LGB_V1_Removable | 3,518,155 (96.1%) | 142,644 (3.9%) | +0.0075 m³/m³ | 0.0484 m³/m³ | 0.0478 m³/m³ | 0.6246 | 0.9667 | +0.0043 [-0.0005, +0.0112] (interprétable) |
| LGB_V2_Rem_ContextRep | 2,848,005 (77.8%) | 812,794 (22.2%) | +0.0059 m³/m³ | 0.0311 m³/m³ | 0.0305 m³/m³ | 0.8638 | 0.8783 | +0.0030 [-0.0033, +0.0106] (interprétable) |
| LGB_V3_Rem_AllContext | 2,787,172 (76.1%) | 873,627 (23.9%) | +0.0049 m³/m³ | 0.0283 m³/m³ | 0.0278 m³/m³ | 0.8910 | 0.8589 | +0.0010 [-0.0028, +0.0062] (interprétable) |

### 5.2 Échantillon commun (1,521,751 observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,521,751 (41.6%) | 0 (0.0%) | -0.0363 m³/m³ | 0.0636 m³/m³ | 0.0522 m³/m³ | 0.5862 | 1.0000 | +0.0107 [-0.0110, +0.0186] (interprétable) |
| B1_Climatology | 1,521,751 (41.6%) | 0 (0.0%) | +0.0113 m³/m³ | 0.0441 m³/m³ | 0.0427 m³/m³ | 0.7461 | 1.0000 | +0.0178 [-0.0094, +0.0271] (interprétable) |
| B2_CDF_Legacy | 1,521,751 (41.6%) | 0 (0.0%) | +0.0087 m³/m³ | 0.0327 m³/m³ | 0.0315 m³/m³ | 0.8237 | 0.8257 | -0.0048 [-0.0135, -0.0000] (interprétable) |
| B2_train_CDF | 1,521,751 (41.6%) | 0 (0.0%) | +0.0134 m³/m³ | 0.0382 m³/m³ | 0.0358 m³/m³ | 0.7766 | 0.8131 | -0.0024 [-0.0118, +0.0018] (interprétable) |
| LGB_V1_Removable | 1,521,751 (41.6%) | 0 (0.0%) | -0.0051 m³/m³ | 0.0486 m³/m³ | 0.0483 m³/m³ | 0.5747 | 0.9700 | +0.0154 [-0.0200, +0.0266] (interprétable) |
| LGB_V2_Rem_ContextRep | 1,521,751 (41.6%) | 0 (0.0%) | +0.0147 m³/m³ | 0.0385 m³/m³ | 0.0356 m³/m³ | 0.8223 | 0.8991 | +0.0036 [-0.0222, +0.0125] (interprétable) |
| LGB_V3_Rem_AllContext | 1,521,751 (41.6%) | 0 (0.0%) | +0.0122 m³/m³ | 0.0348 m³/m³ | 0.0326 m³/m³ | 0.8630 | 0.8801 | -0.0028 [-0.0289, +0.0096] (interprétable) |

### 5.3 Stratification saisonnière des métriques (Échantillon spécifique)

| Référence | Saison | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | DJF | 600,414 (60.4%) | 394,467 (39.6%) | -0.0449 m³/m³ | 0.0684 m³/m³ | 0.0516 m³/m³ | 0.6647 | 1.0000 |
| B0_Raw | MAM | 484,559 (47.1%) | 545,051 (52.9%) | -0.0441 m³/m³ | 0.0645 m³/m³ | 0.0471 m³/m³ | 0.6654 | 1.0000 |
| B0_Raw | JJA | 273,124 (42.8%) | 365,754 (57.2%) | +0.0199 m³/m³ | 0.0388 m³/m³ | 0.0333 m³/m³ | 0.5234 | 1.0000 |
| B0_Raw | SON | 553,975 (55.5%) | 443,455 (44.5%) | -0.0140 m³/m³ | 0.0532 m³/m³ | 0.0513 m³/m³ | 0.6034 | 1.0000 |
| B0_Raw | Annual | 1,912,072 (52.2%) | 1,748,727 (47.8%) | -0.0265 m³/m³ | 0.0597 m³/m³ | 0.0535 m³/m³ | 0.6278 | 1.0000 |
| B1_Climatology | DJF | 991,370 (99.6%) | 3,511 (0.4%) | +0.0143 m³/m³ | 0.0432 m³/m³ | 0.0408 m³/m³ | 0.7466 | 1.0000 |
| B1_Climatology | MAM | 909,843 (88.4%) | 119,767 (11.6%) | +0.0065 m³/m³ | 0.0447 m³/m³ | 0.0442 m³/m³ | 0.7122 | 1.0000 |
| B1_Climatology | JJA | 293,784 (46.0%) | 345,094 (54.0%) | +0.0121 m³/m³ | 0.0282 m³/m³ | 0.0254 m³/m³ | 0.8069 | 1.0000 |
| B1_Climatology | SON | 938,500 (94.1%) | 58,930 (5.9%) | +0.0188 m³/m³ | 0.0439 m³/m³ | 0.0396 m³/m³ | 0.7336 | 1.0000 |
| B1_Climatology | Annual | 3,133,497 (85.6%) | 527,302 (14.4%) | +0.0132 m³/m³ | 0.0427 m³/m³ | 0.0406 m³/m³ | 0.7605 | 1.0000 |
| B2_CDF_Legacy | DJF | 852,868 (85.7%) | 38,523 (3.9%) | +0.0073 m³/m³ | 0.0344 m³/m³ | 0.0336 m³/m³ | 0.7958 | 0.8159 |
| B2_CDF_Legacy | MAM | 793,398 (77.1%) | 83,014 (8.1%) | +0.0047 m³/m³ | 0.0338 m³/m³ | 0.0335 m³/m³ | 0.7894 | 0.8047 |
| B2_CDF_Legacy | JJA | 237,678 (37.2%) | 369,673 (57.9%) | +0.0125 m³/m³ | 0.0302 m³/m³ | 0.0275 m³/m³ | 0.7708 | 0.8114 |
| B2_CDF_Legacy | SON | 738,436 (74.0%) | 161,329 (16.2%) | +0.0130 m³/m³ | 0.0395 m³/m³ | 0.0373 m³/m³ | 0.7383 | 0.7921 |
| B2_CDF_Legacy | Annual | 2,622,380 (71.6%) | 652,539 (17.8%) | +0.0086 m³/m³ | 0.0354 m³/m³ | 0.0343 m³/m³ | 0.8040 | 0.8009 |
| B2_train_CDF | DJF | 891,746 (89.6%) | 20,890 (2.1%) | +0.0172 m³/m³ | 0.0441 m³/m³ | 0.0406 m³/m³ | 0.7184 | 0.8135 |
| B2_train_CDF | MAM | 857,651 (83.3%) | 100,609 (9.8%) | +0.0111 m³/m³ | 0.0413 m³/m³ | 0.0398 m³/m³ | 0.7190 | 0.7682 |
| B2_train_CDF | JJA | 252,249 (39.5%) | 355,653 (55.7%) | +0.0156 m³/m³ | 0.0343 m³/m³ | 0.0305 m³/m³ | 0.7179 | 0.7885 |
| B2_train_CDF | SON | 785,888 (78.8%) | 119,588 (12.0%) | +0.0229 m³/m³ | 0.0469 m³/m³ | 0.0409 m³/m³ | 0.6867 | 0.7775 |
| B2_train_CDF | Annual | 2,787,534 (76.1%) | 596,740 (16.3%) | +0.0168 m³/m³ | 0.0433 m³/m³ | 0.0399 m³/m³ | 0.7451 | 0.7801 |
| LGB_V1_Removable | DJF | 979,670 (98.5%) | 15,211 (1.5%) | -0.0121 m³/m³ | 0.0465 m³/m³ | 0.0449 m³/m³ | 0.6706 | 0.9638 |
| LGB_V1_Removable | MAM | 971,198 (94.3%) | 58,412 (5.7%) | -0.0108 m³/m³ | 0.0474 m³/m³ | 0.0461 m³/m³ | 0.6003 | 0.9714 |
| LGB_V1_Removable | JJA | 595,935 (93.3%) | 42,943 (6.7%) | +0.0453 m³/m³ | 0.0553 m³/m³ | 0.0317 m³/m³ | 0.3555 | 0.9215 |
| LGB_V1_Removable | SON | 971,352 (97.4%) | 26,078 (2.6%) | +0.0223 m³/m³ | 0.0469 m³/m³ | 0.0412 m³/m³ | 0.6750 | 0.9693 |
| LGB_V1_Removable | Annual | 3,518,155 (96.1%) | 142,644 (3.9%) | +0.0075 m³/m³ | 0.0484 m³/m³ | 0.0478 m³/m³ | 0.6246 | 0.9667 |
| LGB_V2_Rem_ContextRep | DJF | 909,912 (91.5%) | 84,969 (8.5%) | +0.0038 m³/m³ | 0.0308 m³/m³ | 0.0306 m³/m³ | 0.8540 | 0.9282 |
| LGB_V2_Rem_ContextRep | MAM | 902,779 (87.7%) | 126,831 (12.3%) | +0.0030 m³/m³ | 0.0307 m³/m³ | 0.0306 m³/m³ | 0.8582 | 0.8644 |
| LGB_V2_Rem_ContextRep | JJA | 261,352 (40.9%) | 377,526 (59.1%) | +0.0048 m³/m³ | 0.0208 m³/m³ | 0.0202 m³/m³ | 0.8634 | 0.7527 |
| LGB_V2_Rem_ContextRep | SON | 773,962 (77.6%) | 223,468 (22.4%) | +0.0123 m³/m³ | 0.0345 m³/m³ | 0.0323 m³/m³ | 0.8501 | 0.8783 |
| LGB_V2_Rem_ContextRep | Annual | 2,848,005 (77.8%) | 812,794 (22.2%) | +0.0059 m³/m³ | 0.0311 m³/m³ | 0.0305 m³/m³ | 0.8638 | 0.8783 |
| LGB_V3_Rem_AllContext | DJF | 905,683 (91.0%) | 89,198 (9.0%) | +0.0030 m³/m³ | 0.0273 m³/m³ | 0.0271 m³/m³ | 0.8900 | 0.9109 |
| LGB_V3_Rem_AllContext | MAM | 903,842 (87.8%) | 125,768 (12.2%) | +0.0047 m³/m³ | 0.0294 m³/m³ | 0.0291 m³/m³ | 0.8873 | 0.8434 |
| LGB_V3_Rem_AllContext | JJA | 236,647 (37.0%) | 402,231 (63.0%) | +0.0026 m³/m³ | 0.0185 m³/m³ | 0.0183 m³/m³ | 0.8854 | 0.7276 |
| LGB_V3_Rem_AllContext | SON | 741,000 (74.3%) | 256,430 (25.7%) | +0.0084 m³/m³ | 0.0305 m³/m³ | 0.0293 m³/m³ | 0.8739 | 0.8627 |
| LGB_V3_Rem_AllContext | Annual | 2,787,172 (76.1%) | 873,627 (23.9%) | +0.0049 m³/m³ | 0.0283 m³/m³ | 0.0278 m³/m³ | 0.8910 | 0.8589 |

---

## 6. Figures de publication (FR et EN, sans interprétation)

- **Importance par groupe (4 groupes) :** [F7_shap_group_importance.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F7_shap_group_importance.png) | EN: [F7_shap_group_importance_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F7_shap_group_importance_en.png)
- **Beeswarm des 15 principales variables :** [F8_shap_beeswarm_top15.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F8_shap_beeswarm_top15.png) | EN: [F8_shap_beeswarm_top15_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F8_shap_beeswarm_top15_en.png)
- **Stabilité des rangs entre plis :** [F9_shap_fold_stability.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F9_shap_fold_stability.png) | EN: [F9_shap_fold_stability_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F9_shap_fold_stability_en.png)
- **Cartes saisonnières de correction moyenne (V1, V2, V3, B2-train) :** [F10_seasonal_applied_corrections_v1_v2_v3.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_seasonal_applied_corrections_v1_v2_v3.png) | EN: [F10_seasonal_applied_corrections_v1_v2_v3_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_seasonal_applied_corrections_v1_v2_v3_en.png)
- **Carte spatiale du masque de rugosité manquante :** [F10_roughness_is_nan.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_roughness_is_nan.png) | EN: [F10_roughness_is_nan_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_roughness_is_nan_en.png)
- **Carte spatiale de la texture dominante Noah-MP :** [F10_dominant_soiltype.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_dominant_soiltype.png) | EN: [F10_dominant_soiltype_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_dominant_soiltype_en.png)
