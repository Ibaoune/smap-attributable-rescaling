# Rapport de Phase 2c : Correction Attribuable LightGBM et TreeSHAP Exact

**Date d'exécution :** 2026-11-23 13:49:07  
**Période d'entraînement :** 2016-01-01 au 2018-12-31 (scindée en 5 plis spatiaux de 1°)  
**Période d'évaluation :** 2019-01-01 au 2019-12-31 (année 2019 scellée au préalable)  
**Année 2020 :** Strictement scellée (aucune lecture, aucune inférence)  
**Population d'entraînement :** Échantillon stratifié de 100 observations par pixel tirées aléatoirement parmi les observations QC-valides (hors flétrissement) sur 2016–2018.  
**Algorithme principal :** LightGBM (300 arbres, learning_rate=0,05, num_leaves=63)  
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
| Pli 0 | 1,132,300 | 292,700 | 0.0213 m³/m³ | 0.9289 | -0.0022 m³/m³ |
| Pli 1 | 1,142,000 | 283,000 | 0.0215 m³/m³ | 0.9328 | -0.0017 m³/m³ |
| Pli 2 | 1,150,800 | 274,200 | 0.0222 m³/m³ | 0.9102 | -0.0006 m³/m³ |
| Pli 3 | 1,135,400 | 289,600 | 0.0202 m³/m³ | 0.9275 | -0.0003 m³/m³ |
| Pli 4 | 1,139,500 | 285,500 | 0.0226 m³/m³ | 0.9086 | +0.0027 m³/m³ |
| **Moyenne** | — | — | **0.0216 m³/m³** | **0.9216** | — |

---

## 3. Stabilité des attributions SHAP au niveau des variables (Jalon H4)

### 3.1 Rangs des 15 variables principales à travers les plis

| Rang | Variable | Rang Moyen | Pli 0 | Pli 1 | Pli 2 | Pli 3 | Pli 4 | SHAP Absolu Moyen |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | `ssm_lag1` | 1.00 | 1 | 1 | 1 | 1 | 1 | 0.0296 |
| 2 | `roughness_coefficient` | 2.20 | 2 | 2 | 2 | 2 | 3 | 0.0092 |
| 3 | `vegetation_water_content` | 3.60 | 3 | 3 | 3 | 4 | 5 | 0.0062 |
| 4 | `doy_sin` | 4.20 | 4 | 4 | 4 | 5 | 4 | 0.0058 |
| 5 | `ssm_rel_wilt_ref_lag1` | 4.20 | 6 | 5 | 5 | 3 | 2 | 0.0062 |
| 6 | `pass_is_pm` | 6.20 | 7 | 6 | 6 | 6 | 6 | 0.0044 |
| 7 | `vegetation_opacity` | 8.20 | 8 | 8 | 9 | 7 | 9 | 0.0037 |
| 8 | `surface_temperature` | 8.40 | 10 | 9 | 8 | 8 | 7 | 0.0036 |
| 9 | `roughness_is_nan` | 9.40 | 5 | 7 | 11 | 9 | 15 | 0.0035 |
| 10 | `lai_lag1` | 10.00 | 9 | 10 | 10 | 11 | 10 | 0.0034 |
| 11 | `sm_layer4_lag1` | 10.20 | 12 | 11 | 7 | 10 | 11 | 0.0030 |
| 12 | `sm_layer3_lag1` | 11.40 | 11 | 12 | 13 | 13 | 8 | 0.0028 |
| 13 | `et_lag1` | 13.80 | 14 | 15 | 12 | 12 | 16 | 0.0023 |
| 14 | `precip_mswep` | 13.80 | 13 | 13 | 15 | 14 | 14 | 0.0021 |
| 15 | `doy_cos` | 15.00 | 15 | 14 | 14 | 15 | 17 | 0.0020 |

### 3.2 Matrice de corrélation de Spearman entre plis

| Paire de plis | Spearman $ho$ (Top 15) | Spearman $ho$ (Toutes les 66 variables) |
| :--- | :--- | :--- |
| Pli 0 vs Pli 1 | 0.9786 (p=2.65e-10) | 0.9692 (p=1.12e-40) |
| Pli 0 vs Pli 2 | 0.8536 (p=5.17e-05) | 0.9819 (p=5.91e-48) |
| Pli 0 vs Pli 3 | 0.9107 (p=2.40e-06) | 0.9871 (p=1.35e-52) |
| Pli 0 vs Pli 4 | 0.8071 (p=2.75e-04) | 0.9814 (p=1.40e-47) |
| Pli 1 vs Pli 2 | 0.9143 (p=1.85e-06) | 0.9713 (p=1.27e-41) |
| Pli 1 vs Pli 3 | 0.9536 (p=3.79e-08) | 0.9733 (p=1.24e-42) |
| Pli 1 vs Pli 4 | 0.8679 (p=2.75e-05) | 0.9583 (p=1.64e-36) |
| Pli 2 vs Pli 3 | 0.9536 (p=3.79e-08) | 0.9876 (p=3.10e-53) |
| Pli 2 vs Pli 4 | 0.8679 (p=2.75e-05) | 0.9686 (p=2.20e-40) |
| Pli 3 vs Pli 4 | 0.8929 (p=7.49e-06) | 0.9761 (p=3.93e-44) |
| **Moyenne** | **0.9000** | **0.9755** |

---

## 4. Variantes de correction et décomposition saisonnière (2019)

Les 3 variantes de correction appliquent la décomposition additive :
- **V1 (removable) :** $\mathrm{SMAP}^*_{\mathrm{V1}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \phi_0)$
- **V2 (removable + context_rep) :** $\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$
- **V3 (removable + all context) :** $\mathrm{SMAP}^*_{\mathrm{V3}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \Phi_{\mathrm{context\_proc}} + \phi_0)$

### Correction moyenne appliquée $\overline{\Delta \mathrm{SMAP}}$ et biais résiduel $\overline{d_{\mathrm{res}}}$ par saison

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{\mathrm{res}}$ | V2 : Correction | V2 : Biais $d_{\mathrm{res}}$ | V3 : Correction | V3 : Biais $d_{\mathrm{res}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| DJF | 994,881 | +0.0400 m³/m³ | -0.0137 m³/m³ | +0.0564 m³/m³ | +0.0027 m³/m³ | +0.0558 m³/m³ | +0.0021 m³/m³ |
| MAM | 1,029,610 | +0.0448 m³/m³ | -0.0140 m³/m³ | +0.0607 m³/m³ | +0.0019 m³/m³ | +0.0622 m³/m³ | +0.0034 m³/m³ |
| JJA | 638,878 | +0.0454 m³/m³ | +0.0402 m³/m³ | +0.0052 m³/m³ | -0.0000 m³/m³ | +0.0029 m³/m³ | -0.0023 m³/m³ |
| SON | 997,430 | +0.0469 m³/m³ | +0.0206 m³/m³ | +0.0350 m³/m³ | +0.0087 m³/m³ | +0.0313 m³/m³ | +0.0050 m³/m³ |
| Annual | 3,660,799 | +0.0442 m³/m³ | +0.0050 m³/m³ | +0.0428 m³/m³ | +0.0036 m³/m³ | +0.0417 m³/m³ | +0.0025 m³/m³ |

---

## 5. Protocole 2a complet : Comparaison avec B0, B1, B2, B2-train et RF

Le flétrissement ($\mathrm{SMAP}^* \ge 	heta_{\mathrm{wlt}} + 0{,}02$) est recalculé sur l'observation corrigée pour chaque référence.

### 5.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | Pearson $r$ | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,912,072 (52.2%) | 1,748,727 (47.8%) | -0.0265 m³/m³ | 0.0597 m³/m³ | 0.0535 m³/m³ | 0.6278 | +0.0109 m³/m³ |
| B1_Climatology | 3,096,112 (84.6%) | 564,687 (15.4%) | +0.0144 m³/m³ | 0.0446 m³/m³ | 0.0423 m³/m³ | 0.7440 | +0.0049 m³/m³ |
| B2_CDF_Legacy | 2,622,380 (71.6%) | 652,539 (17.8%) | +0.0086 m³/m³ | 0.0354 m³/m³ | 0.0343 m³/m³ | 0.8040 | +0.0010 m³/m³ |
| B2_train_CDF | 2,787,534 (76.1%) | 596,740 (16.3%) | +0.0168 m³/m³ | 0.0433 m³/m³ | 0.0399 m³/m³ | 0.7451 | +0.0026 m³/m³ |
| LGB_V1_Removable | 3,518,886 (96.1%) | 141,913 (3.9%) | +0.0073 m³/m³ | 0.0485 m³/m³ | 0.0479 m³/m³ | 0.6227 | +0.0037 m³/m³ |
| LGB_V2_Rem_ContextRep | 2,843,171 (77.7%) | 817,628 (22.3%) | +0.0057 m³/m³ | 0.0310 m³/m³ | 0.0305 m³/m³ | 0.8643 | +0.0019 m³/m³ |
| LGB_V3_Rem_AllContext | 2,787,556 (76.1%) | 873,243 (23.9%) | +0.0049 m³/m³ | 0.0284 m³/m³ | 0.0279 m³/m³ | 0.8902 | +0.0007 m³/m³ |
| RF_Comparison_Phase2b | 3,548,975 (96.9%) | 111,824 (3.1%) | +0.0061 m³/m³ | 0.0456 m³/m³ | 0.0452 m³/m³ | 0.7729 | N/A m³/m³ |

### 5.2 Échantillon commun (1,515,650 observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | Pearson $r$ | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,515,650 (41.4%) | 0 (0.0%) | -0.0364 m³/m³ | 0.0636 m³/m³ | 0.0522 m³/m³ | 0.5848 | +0.0108 m³/m³ |
| B1_Climatology | 1,515,650 (41.4%) | 0 (0.0%) | +0.0119 m³/m³ | 0.0460 m³/m³ | 0.0444 m³/m³ | 0.7296 | +0.0167 m³/m³ |
| B2_CDF_Legacy | 1,515,650 (41.4%) | 0 (0.0%) | +0.0088 m³/m³ | 0.0328 m³/m³ | 0.0316 m³/m³ | 0.8223 | -0.0050 m³/m³ |
| B2_train_CDF | 1,515,650 (41.4%) | 0 (0.0%) | +0.0135 m³/m³ | 0.0383 m³/m³ | 0.0358 m³/m³ | 0.7749 | -0.0025 m³/m³ |
| LGB_V1_Removable | 1,515,650 (41.4%) | 0 (0.0%) | -0.0056 m³/m³ | 0.0486 m³/m³ | 0.0483 m³/m³ | 0.5744 | +0.0124 m³/m³ |
| LGB_V2_Rem_ContextRep | 1,515,650 (41.4%) | 0 (0.0%) | +0.0145 m³/m³ | 0.0383 m³/m³ | 0.0355 m³/m³ | 0.8226 | +0.0004 m³/m³ |
| LGB_V3_Rem_AllContext | 1,515,650 (41.4%) | 0 (0.0%) | +0.0124 m³/m³ | 0.0350 m³/m³ | 0.0327 m³/m³ | 0.8614 | -0.0064 m³/m³ |

---

## 6. Figures de publication (FR et EN, sans interprétation)

- **Importance par groupe (4 groupes) :** [F7_shap_group_importance.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F7_shap_group_importance.png) | EN: [F7_shap_group_importance_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F7_shap_group_importance_en.png)
- **Beeswarm des 15 principales variables :** [F8_shap_beeswarm_top15.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F8_shap_beeswarm_top15.png) | EN: [F8_shap_beeswarm_top15_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F8_shap_beeswarm_top15_en.png)
- **Stabilité des rangs entre plis :** [F9_shap_fold_stability.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F9_shap_fold_stability.png) | EN: [F9_shap_fold_stability_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F9_shap_fold_stability_en.png)
- **Cartes saisonnières de correction moyenne (V1, V2, V3) :** [F10_seasonal_applied_corrections_v1_v2_v3.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F10_seasonal_applied_corrections_v1_v2_v3.png) | EN: [F10_seasonal_applied_corrections_v1_v2_v3_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F10_seasonal_applied_corrections_v1_v2_v3_en.png)
