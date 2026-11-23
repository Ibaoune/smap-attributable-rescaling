# Rapport de Phase 2c : Correction Attribuable LightGBM et TreeSHAP Exact

**Date d'exécution :** 2026-11-23 16:27:59  
**Période d'entraînement :** 2016-01-01 au 2018-12-31 (scindée en 5 plis spatiaux de 1°)  
**Période d'évaluation :** 2019-01-01 au 2019-12-31 (année 2019 scellée au préalable)  
**Année 2020 :** Strictement scellée (aucune lecture, aucune inférence)  
**Population d'entraînement :** Échantillon stratifié de 100 observations par pixel tirées aléatoirement parmi les observations QC-valides (hors flétrissement) sur 2016–2018.  
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

Validation sur blocs spatiaux disjoints de 1° avec calcul du TreeSHAP exact sur au moins 25 000 observations de validation par pli :

| Pli | Taille Train | Taille Val | RMSE Val | Pearson r | Biais Val |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Pli 0 | 1,152,800 | 272,200 | 0.0178 m³/m³ | 0.9514 | +0.0037 m³/m³ |
| Pli 1 | 1,095,100 | 329,900 | 0.0230 m³/m³ | 0.9181 | -0.0021 m³/m³ |
| Pli 2 | 1,150,900 | 274,100 | 0.0249 m³/m³ | 0.8885 | -0.0006 m³/m³ |
| Pli 3 | 1,154,400 | 270,600 | 0.0227 m³/m³ | 0.9254 | -0.0045 m³/m³ |
| Pli 4 | 1,146,800 | 278,200 | 0.0181 m³/m³ | 0.9420 | +0.0027 m³/m³ |
| **Moyenne** | — | — | **0.0213 m³/m³** | **0.9251** | — |

---

## 3. Stabilité des attributions SHAP (Groupes, Intra-groupe et Variables)

### 3.1 Stabilité au niveau des 4 GROUPES et corrélation de rang intra-groupe

| Groupe | Contribution Absolue Moyenne | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | Rang Moyen | Spearman intra-groupe moyen |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `removable` | 0.0190 m³/m³ | 2 | 2 | 2 | 2 | 2 | 2.00 | 0.9722 |
| `retained` | 0.0054 m³/m³ | 4 | 4 | 4 | 4 | 4 | 4.00 | 1.0000 |
| `context_representation` | 0.0378 m³/m³ | 1 | 1 | 1 | 1 | 1 | 1.00 | 0.8800 |
| `context_process` | 0.0075 m³/m³ | 3 | 3 | 3 | 3 | 3 | 3.00 | 0.8452 |
| **Moyenne groupes** | — | — | — | — | — | — | — | **1.0000** |

### 3.2 Rangs des 15 variables principales à travers les plis

| Rang | Variable | Rang Moyen | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | SHAP Absolu Moyen |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | `ssm_lag1` | 1.00 | 1 | 1 | 1 | 1 | 1 | 0.0314 |
| 2 | `roughness_coefficient` | 2.00 | 2 | 2 | 2 | 2 | 2 | 0.0127 |
| 3 | `vegetation_water_content` | 3.60 | 3 | 4 | 3 | 3 | 5 | 0.0059 |
| 4 | `doy_sin` | 4.00 | 5 | 3 | 4 | 4 | 4 | 0.0055 |
| 5 | `ssm_rel_wilt_ref_lag1` | 5.20 | 4 | 9 | 5 | 5 | 3 | 0.0056 |
| 6 | `pass_is_pm` | 5.80 | 6 | 5 | 6 | 6 | 6 | 0.0045 |
| 7 | `vegetation_opacity` | 7.00 | 7 | 6 | 7 | 8 | 7 | 0.0037 |
| 8 | `surface_temperature` | 8.00 | 10 | 7 | 8 | 7 | 8 | 0.0033 |
| 9 | `sm_layer3_lag1` | 9.40 | 8 | 11 | 9 | 10 | 9 | 0.0030 |
| 10 | `lai_lag1` | 10.20 | 9 | 8 | 15 | 9 | 10 | 0.0030 |
| 11 | `sm_layer4_lag1` | 10.80 | 12 | 10 | 10 | 11 | 11 | 0.0027 |
| 12 | `et_lag1` | 12.00 | 11 | 12 | 11 | 14 | 12 | 0.0023 |
| 13 | `precip_mswep` | 13.00 | 13 | 14 | 12 | 13 | 13 | 0.0020 |
| 14 | `doy_cos` | 14.60 | 14 | 15 | 13 | 15 | 16 | 0.0019 |
| 15 | `precip_3d_sum` | 15.80 | 15 | 17 | 16 | 16 | 15 | 0.0017 |

### 3.3 Matrice de corrélation de Spearman entre plis

| Paire de plis | Spearman $\rho$ (4 Groupes) | Spearman $\rho$ (Top 15) | Spearman $\rho$ (Toutes les 66 variables) |
| :--- | :--- | :--- | :--- |
| Pli 0 vs Pli 1 | 1.0000 | 0.9000 (p=4.87e-06) | 0.9819 (p=5.50e-48) |
| Pli 0 vs Pli 2 | 1.0000 | 0.9321 (p=4.24e-07) | 0.9843 (p=6.22e-50) |
| Pli 0 vs Pli 3 | 1.0000 | 0.9607 (p=1.30e-08) | 0.9860 (p=1.75e-51) |
| Pli 0 vs Pli 4 | 1.0000 | 0.9714 (p=1.69e-09) | 0.9820 (p=5.11e-48) |
| Pli 1 vs Pli 2 | 1.0000 | 0.8857 (p=1.12e-05) | 0.9838 (p=1.83e-49) |
| Pli 1 vs Pli 3 | 1.0000 | 0.9500 (p=6.09e-08) | 0.9853 (p=8.38e-51) |
| Pli 1 vs Pli 4 | 1.0000 | 0.9071 (p=3.06e-06) | 0.9803 (p=7.92e-47) |
| Pli 2 vs Pli 3 | 1.0000 | 0.9393 (p=2.09e-07) | 0.9790 (p=5.97e-46) |
| Pli 2 vs Pli 4 | 1.0000 | 0.9429 (p=1.42e-07) | 0.9699 (p=5.34e-41) |
| Pli 3 vs Pli 4 | 1.0000 | 0.9714 (p=1.69e-09) | 0.9818 (p=6.84e-48) |
| **Moyenne** | **1.0000** | **0.9361** | **0.9814** |

---

## 4. Variantes de correction et décomposition saisonnière (2019)

Les 3 variantes de correction appliquent la décomposition additive :
- **V1 (removable) :** $\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$
- **V2 (removable + context_rep) :** $\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$
- **V3 (removable + all context) :** $\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$

### Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ et biais résiduel $\overline{d_{\mathrm{res}}}$ par saison

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{\mathrm{res}}$ | V2 : Correction | V2 : Biais $d_{\mathrm{res}}$ | V3 : Correction | V3 : Biais $d_{\mathrm{res}}$ | B2-train : Correction | B2-train : Biais $d_{\mathrm{res}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| DJF | 994,881 | +0.0405 m³/m³ | -0.0132 m³/m³ | +0.0568 m³/m³ | +0.0031 m³/m³ | +0.0558 m³/m³ | +0.0022 m³/m³ | +0.0709 m³/m³ | +0.0162 m³/m³ |
| MAM | 1,029,610 | +0.0448 m³/m³ | -0.0140 m³/m³ | +0.0608 m³/m³ | +0.0020 m³/m³ | +0.0622 m³/m³ | +0.0034 m³/m³ | +0.0652 m³/m³ | +0.0046 m³/m³ |
| JJA | 638,878 | +0.0451 m³/m³ | +0.0399 m³/m³ | +0.0052 m³/m³ | +0.0000 m³/m³ | +0.0028 m³/m³ | -0.0024 m³/m³ | +0.0061 m³/m³ | +0.0007 m³/m³ |
| SON | 997,430 | +0.0469 m³/m³ | +0.0206 m³/m³ | +0.0352 m³/m³ | +0.0089 m³/m³ | +0.0312 m³/m³ | +0.0049 m³/m³ | +0.0441 m³/m³ | +0.0170 m³/m³ |
| Annual | 3,660,799 | +0.0443 m³/m³ | +0.0051 m³/m³ | +0.0430 m³/m³ | +0.0038 m³/m³ | +0.0416 m³/m³ | +0.0025 m³/m³ | +0.0505 m³/m³ | +0.0103 m³/m³ |

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
| LGB_V1_Removable | 3,517,036 (96.1%) | 143,763 (3.9%) | +0.0074 m³/m³ | 0.0484 m³/m³ | 0.0478 m³/m³ | 0.6248 | 0.9667 | +0.0041 [-0.0005, +0.0108] (interprétable) |
| LGB_V2_Rem_ContextRep | 2,843,830 (77.7%) | 816,969 (22.3%) | +0.0059 m³/m³ | 0.0310 m³/m³ | 0.0305 m³/m³ | 0.8640 | 0.8780 | +0.0029 [-0.0030, +0.0105] (interprétable) |
| LGB_V3_Rem_AllContext | 2,787,716 (76.2%) | 873,083 (23.8%) | +0.0049 m³/m³ | 0.0283 m³/m³ | 0.0279 m³/m³ | 0.8903 | 0.8596 | +0.0019 [-0.0019, +0.0068] (interprétable) |

### 5.2 Échantillon commun (1,523,435 observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,523,435 (41.6%) | 0 (0.0%) | -0.0362 m³/m³ | 0.0635 m³/m³ | 0.0522 m³/m³ | 0.5863 | 1.0000 | +0.0109 [-0.0113, +0.0187] (interprétable) |
| B1_Climatology | 1,523,435 (41.6%) | 0 (0.0%) | +0.0113 m³/m³ | 0.0441 m³/m³ | 0.0427 m³/m³ | 0.7465 | 1.0000 | +0.0177 [-0.0101, +0.0270] (interprétable) |
| B2_CDF_Legacy | 1,523,435 (41.6%) | 0 (0.0%) | +0.0087 m³/m³ | 0.0327 m³/m³ | 0.0316 m³/m³ | 0.8238 | 0.8254 | -0.0049 [-0.0144, -0.0001] (interprétable) |
| B2_train_CDF | 1,523,435 (41.6%) | 0 (0.0%) | +0.0134 m³/m³ | 0.0383 m³/m³ | 0.0358 m³/m³ | 0.7767 | 0.8129 | -0.0024 [-0.0128, +0.0017] (interprétable) |
| LGB_V1_Removable | 1,523,435 (41.6%) | 0 (0.0%) | -0.0053 m³/m³ | 0.0486 m³/m³ | 0.0483 m³/m³ | 0.5749 | 0.9701 | +0.0150 [-0.0219, +0.0269] (interprétable) |
| LGB_V2_Rem_ContextRep | 1,523,435 (41.6%) | 0 (0.0%) | +0.0146 m³/m³ | 0.0384 m³/m³ | 0.0355 m³/m³ | 0.8235 | 0.8990 | +0.0026 [-0.0251, +0.0126] (interprétable) |
| LGB_V3_Rem_AllContext | 1,523,435 (41.6%) | 0 (0.0%) | +0.0121 m³/m³ | 0.0348 m³/m³ | 0.0327 m³/m³ | 0.8622 | 0.8812 | -0.0036 [-0.0308, +0.0092] (interprétable) |

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
| LGB_V1_Removable | DJF | 979,579 (98.5%) | 15,302 (1.5%) | -0.0123 m³/m³ | 0.0465 m³/m³ | 0.0449 m³/m³ | 0.6698 | 0.9635 |
| LGB_V1_Removable | MAM | 970,636 (94.3%) | 58,974 (5.7%) | -0.0108 m³/m³ | 0.0474 m³/m³ | 0.0461 m³/m³ | 0.6007 | 0.9714 |
| LGB_V1_Removable | JJA | 595,862 (93.3%) | 43,016 (6.7%) | +0.0453 m³/m³ | 0.0552 m³/m³ | 0.0316 m³/m³ | 0.3605 | 0.9240 |
| LGB_V1_Removable | SON | 970,959 (97.3%) | 26,471 (2.7%) | +0.0223 m³/m³ | 0.0468 m³/m³ | 0.0411 m³/m³ | 0.6766 | 0.9694 |
| LGB_V1_Removable | Annual | 3,517,036 (96.1%) | 143,763 (3.9%) | +0.0074 m³/m³ | 0.0484 m³/m³ | 0.0478 m³/m³ | 0.6248 | 0.9667 |
| LGB_V2_Rem_ContextRep | DJF | 910,148 (91.5%) | 84,733 (8.5%) | +0.0036 m³/m³ | 0.0308 m³/m³ | 0.0305 m³/m³ | 0.8538 | 0.9283 |
| LGB_V2_Rem_ContextRep | MAM | 901,947 (87.6%) | 127,663 (12.4%) | +0.0032 m³/m³ | 0.0309 m³/m³ | 0.0307 m³/m³ | 0.8579 | 0.8645 |
| LGB_V2_Rem_ContextRep | JJA | 257,114 (40.2%) | 381,764 (59.8%) | +0.0046 m³/m³ | 0.0207 m³/m³ | 0.0201 m³/m³ | 0.8623 | 0.7501 |
| LGB_V2_Rem_ContextRep | SON | 774,621 (77.7%) | 222,809 (22.3%) | +0.0122 m³/m³ | 0.0342 m³/m³ | 0.0320 m³/m³ | 0.8522 | 0.8775 |
| LGB_V2_Rem_ContextRep | Annual | 2,843,830 (77.7%) | 816,969 (22.3%) | +0.0059 m³/m³ | 0.0310 m³/m³ | 0.0305 m³/m³ | 0.8640 | 0.8780 |
| LGB_V3_Rem_AllContext | DJF | 905,247 (91.0%) | 89,634 (9.0%) | +0.0028 m³/m³ | 0.0273 m³/m³ | 0.0272 m³/m³ | 0.8887 | 0.9112 |
| LGB_V3_Rem_AllContext | MAM | 902,812 (87.7%) | 126,798 (12.3%) | +0.0047 m³/m³ | 0.0296 m³/m³ | 0.0292 m³/m³ | 0.8857 | 0.8451 |
| LGB_V3_Rem_AllContext | JJA | 235,661 (36.9%) | 403,217 (63.1%) | +0.0028 m³/m³ | 0.0184 m³/m³ | 0.0182 m³/m³ | 0.8857 | 0.7269 |
| LGB_V3_Rem_AllContext | SON | 743,996 (74.6%) | 253,434 (25.4%) | +0.0084 m³/m³ | 0.0304 m³/m³ | 0.0292 m³/m³ | 0.8751 | 0.8624 |
| LGB_V3_Rem_AllContext | Annual | 2,787,716 (76.2%) | 873,083 (23.8%) | +0.0049 m³/m³ | 0.0283 m³/m³ | 0.0279 m³/m³ | 0.8903 | 0.8596 |

---

## 6. Figures de publication (FR et EN, sans interprétation)

- **Importance par groupe (4 groupes) :** [F7_shap_group_importance.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F7_shap_group_importance.png) | EN: [F7_shap_group_importance_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F7_shap_group_importance_en.png)
- **Beeswarm des 15 principales variables :** [F8_shap_beeswarm_top15.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F8_shap_beeswarm_top15.png) | EN: [F8_shap_beeswarm_top15_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F8_shap_beeswarm_top15_en.png)
- **Stabilité des rangs entre plis :** [F9_shap_fold_stability.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F9_shap_fold_stability.png) | EN: [F9_shap_fold_stability_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F9_shap_fold_stability_en.png)
- **Cartes saisonnières de correction moyenne (V1, V2, V3, B2-train) :** [F10_seasonal_applied_corrections_v1_v2_v3.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_seasonal_applied_corrections_v1_v2_v3.png) | EN: [F10_seasonal_applied_corrections_v1_v2_v3_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_seasonal_applied_corrections_v1_v2_v3_en.png)
- **Carte spatiale du masque de rugosité manquante :** [F10_roughness_is_nan.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_roughness_is_nan.png) | EN: [F10_roughness_is_nan_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_roughness_is_nan_en.png)
- **Carte spatiale de la texture dominante Noah-MP :** [F10_dominant_soiltype.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_dominant_soiltype.png) | EN: [F10_dominant_soiltype_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_dominant_soiltype_en.png)
