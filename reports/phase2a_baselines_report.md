# Rapport de Synthèse Phase 2a : Évaluation des Références (Baselines B0, B1, B2)

**Date d'exécution** : 2026-11-23 11:42:22  
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
   - Assimilation de l'observation brute : $\text{smap}^*_0 = \text{smap_sm}$, avec $d_{\text{res}, 0} = \text{smap_sm} - \text{ssm_ol_daily}$.
2. **B1 (Climatologie pixel-mois)** :
   - Moyenne de $d$ par pixel et par mois, $\hat{c}(p, m)$, calculée sur la période 2016–2018 pour toutes les observations passant les règles de contrôle qualité hors flétrissement.
   - Sur 2019 : $\text{smap}^*_1 = \text{smap_sm} - \hat{c}(p, m)$ et $d_{\text{res}, 1} = d - \hat{c}(p, m)$.
3. **B2 (CDF matching v2 hors ligne)** :
   - Application de la transformation `LIS_rescale_with_CDF_matching` avec les fichiers LDT `noahmp_model_cdf_NorthMor.nc` et `smap_obs_cdf_NorthMor.nc`.
   - Sur 2019 : $\text{smap}^*_2 = \text{CDF_rescale}(\text{smap_sm}, p, m)$ et $d_{\text{res}, 2} = \text{smap}^*_2 - \text{ssm_ol_daily}$.

---

## 2. Vérification Numérique du CDF Matching contre DAOBS LIS (Record 2)

Confrontation de la fonction hors ligne avec l'enregistrement 2 des DAOBS sur 5 dates échantillonnées :

| Date et Heure DAOBS | Nombre de points comparés | Biais absolu moyen | Différence absolue max | Concordance exacte (< 1e-5) |
|:---|:---:|:---:|:---:|:---:|
| `2016-01-01 06:00` | 5,420 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-01-01 06:00` | 4,011 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-04-01 06:00` | 8,532 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-07-24 06:00` | 398 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-10-01 18:00` | 647 | 0.00000000 | 0.00000000 | 100.00% |

---

## 3. Données Non Mises à l'Échelle par B2 (IQR et Hors Plage)

Nombre et proportion d'observations du socle hors flétrissement (3,660,799 observations) pour lesquelles la transformation CDF B2 ne renvoie pas de valeur numérique valide (critère IQR < 0.05 de la plage dynamique, ou valeur en dehors de l'intervalle [0.01, 1.0]) :

| Période / Saison | Observations socle (hors flétrissement) | B2 valeurs valides | B2 part valide (%) | B2 non transformé (NaN / IQR) | B2 part non transformée (%) |
|:---|:---:|:---:|:---:|:---:|:---:|
| All | 3,660,799 | 3,274,919 | 89.46% | 385,880 | 10.54% |
| DJF | 994,881 | 891,391 | 89.60% | 103,490 | 10.40% |
| MAM | 1,029,610 | 876,412 | 85.12% | 153,198 | 14.88% |
| JJA | 638,878 | 607,351 | 95.07% | 31,527 | 4.93% |
| SON | 997,430 | 899,765 | 90.21% | 97,665 | 9.79% |

---

## 4. Population d'Évaluation et Observations Assimilables Après Correction

Pour chaque référence, le seuil de flétrissement ($\text{smap}^* < \theta_{\text{wlt}} + 0.02$) est recalculé sur l'observation issue de la correction. Une observation est assimilable si elle satisfait les contrôles hors flétrissement, possède une valeur corrigée valide (non-NaN, comprise entre 0.01 et 1.0) et respecte $\text{smap}^* \ge \theta_{\text{wlt}} + 0.02$ :

| Échantillon | Référence | Observations assimilables | Part du socle hors flétrissement (3,660,799 obs) | Part du total 2019 (3,815,354 obs) |
|:---|:---|:---:|:---:|:---:|
| Spécifique | **B0 (NoCDF)** | 1,912,072 | 52.23% | 50.12% |
| Spécifique | **B1 (Climatologie)** | 3,133,887 | 85.61% | 82.14% |
| Spécifique | **B2 (CDF v2)** | 2,622,380 | 71.63% | 68.73% |
| Commun | **B0 / B1 / B2** | 1,615,091 | 44.12% | 42.33% |

---

## 5. Performances Globales sur l'Année 2019

Calculées sur l'échantillon assimilable propre à chaque méthode, ainsi que sur l'échantillon commun aux trois références (1,615,091 observations) :

| Référence | Biais résiduel (Échantillon spécifique) | RMSE résiduelle (Échantillon spécifique) | Biais résiduel (Échantillon commun) | RMSE résiduelle (Échantillon commun) |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | -0.0265 m³/m³ | 0.0597 m³/m³ | -0.0329 m³/m³ | 0.0623 m³/m³ |
| **B1 (Climatologie)** | +0.0135 m³/m³ | 0.0427 m³/m³ | +0.0127 m³/m³ | 0.0440 m³/m³ |
| **B2 (CDF v2)** | +0.0086 m³/m³ | 0.0354 m³/m³ | +0.0100 m³/m³ | 0.0336 m³/m³ |

- **B0** : Biais de -0.0265 m³/m³ (spécifique) et -0.0329 m³/m³ (commun) ; RMSE de 0.0597 m³/m³ (spécifique) et 0.0623 m³/m³ (commun).
- **B1** : Biais de +0.0135 m³/m³ (spécifique) et +0.0127 m³/m³ (commun) ; RMSE de 0.0427 m³/m³ (spécifique) et 0.0440 m³/m³ (commun).
- **B2** : Biais de +0.0086 m³/m³ (spécifique) et +0.0100 m³/m³ (commun) ; RMSE de 0.0354 m³/m³ (spécifique) et 0.0336 m³/m³ (commun).

---

## 6. Corrélation d'Anomalies avec SMAP Brut

Corrélation de Pearson entre les anomalies temporelles de l'observation corrigée et celles de SMAP brut (anomalie = observation moins moyenne locale du pixel pour le mois considéré) :

| Référence | Anomaly Correlation (Échantillon spécifique) | N observations spécifique | Anomaly Correlation (Échantillon commun) | N observations commun |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | 1.0000 | 1,912,072 | 1.0000 | 1,615,091 |
| **B1 (Climatologie)** | 1.0000 | 3,133,887 | 1.0000 | 1,615,091 |
| **B2 (CDF v2)** | 0.8009 | 2,622,380 | 0.8173 | 1,615,091 |

---

## 7. Résultats Stratifiés par Saison (Année 2019)

| Saison | RMSE B0 (spécifique) | RMSE B1 (spécifique) | RMSE B2 (spécifique) | N Obs Commun | RMSE B0 (commun) | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| DJF | 0.0684 | 0.0433 | 0.0344 | 569,762 | 0.0692 | 0.0409 | 0.0313 |
| MAM | 0.0645 | 0.0447 | 0.0338 | 433,639 | 0.0661 | 0.0495 | 0.0318 |
| JJA | 0.0388 | 0.0282 | 0.0302 | 132,609 | 0.0397 | 0.0332 | 0.0345 |
| SON | 0.0532 | 0.0439 | 0.0395 | 479,081 | 0.0549 | 0.0450 | 0.0375 |

---

## 8. Résultats Stratifiés par Texture Dominante (Top 5 Classes)

| Texture | N Obs B0 | RMSE B0 | N Obs B1 | RMSE B1 | N Obs B2 | RMSE B2 | N Obs Commun | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Loam | 1,430,068 | 0.0613 | 2,364,879 | 0.0441 | 1,919,404 | 0.0371 | 1,223,936 | 0.0439 | 0.0342 |
| Clay Loam | 98,285 | 0.0702 | 252,273 | 0.0378 | 248,232 | 0.0317 | 96,962 | 0.0464 | 0.0350 |
| Sandy Clay Loam | 137,156 | 0.0449 | 152,784 | 0.0413 | 140,089 | 0.0291 | 120,083 | 0.0441 | 0.0293 |
| Sandy Loam | 53,757 | 0.0699 | 102,727 | 0.0427 | 96,303 | 0.0326 | 51,622 | 0.0524 | 0.0343 |
| Clay | 46,540 | 0.0658 | 121,089 | 0.0366 | 111,842 | 0.0251 | 41,721 | 0.0500 | 0.0333 |

---

## 9. Résultats Stratifiés par Classe d'Irrigation GMIA

| Classe GMIA | N Obs B0 | RMSE B0 | N Obs B1 | RMSE B1 | N Obs B2 | RMSE B2 | N Obs Commun | RMSE B1 (commun) | RMSE B2 (commun) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0% | 1,211,241 | 0.0578 | 2,070,794 | 0.0439 | 1,668,791 | 0.0366 | 1,006,440 | 0.0446 | 0.0341 |
| 0-10% | 601,277 | 0.0622 | 849,679 | 0.0405 | 763,274 | 0.0324 | 522,522 | 0.0427 | 0.0321 |
| 10-30% | 56,523 | 0.0667 | 113,588 | 0.0410 | 97,233 | 0.0362 | 46,753 | 0.0448 | 0.0370 |
| >30% | 43,031 | 0.0652 | 99,826 | 0.0383 | 93,082 | 0.0354 | 39,376 | 0.0448 | 0.0380 |

---

## 10. Contraste du Résidu entre GMIA >30 % et GMIA 0 % à Texture Fixée

Écart $\Delta d_{\text{res}} = \bar{d}_{\text{res}}(\text{GMIA} > 30\%) - \bar{d}_{\text{res}}(\text{GMIA} = 0\%)$ calculé à texture fixée (Clay Loam et Loam), avec intervalle de confiance à 95 % par rééchantillonnage bootstrap sur blocs spatiaux de 1° ($B=500$) et nombre de pixels distincts :

| Texture | Saison | Échantillon | Contraste B0 [IC 95%] (px >30%/0%) | Contraste B1 [IC 95%] (px >30%/0%) | Contraste B2 [IC 95%] (px >30%/0%) |
|:---|:---|:---|:---:|:---:|:---:|
| Clay Loam | All | specific | -0.0025 [-0.0169, +0.0103] (252/331 px, 9/15 blk) | -0.0049 [-0.0076, -0.0019] (252/331 px, 9/15 blk) | -0.0063 [-0.0123, -0.0027] (249/328 px, 9/15 blk) |
| Clay Loam | All | common | -0.0034 [-0.0192, +0.0090] (249/328 px, 9/15 blk) | -0.0128 [-0.0212, -0.0042] (249/328 px, 9/15 blk) | -0.0094 [-0.0190, +0.0018] (249/328 px, 9/15 blk) |
| Clay Loam | DJF | specific | -0.0010 [-0.0198, +0.0150] (252/331 px, 9/15 blk) | -0.0105 [-0.0204, -0.0028] (252/331 px, 9/15 blk) | -0.0074 [-0.0166, -0.0003] (249/328 px, 9/15 blk) |
| Clay Loam | DJF | common | -0.0020 [-0.0224, +0.0152] (249/328 px, 9/15 blk) | -0.0164 [-0.0306, -0.0043] (249/328 px, 9/15 blk) | -0.0105 [-0.0244, +0.0036] (249/328 px, 9/15 blk) |
| Clay Loam | MAM | specific | -0.0062 [-0.0146, +0.0014] (252/331 px, 9/15 blk) | -0.0032 [-0.0113, +0.0048] (252/331 px, 9/15 blk) | -0.0096 [-0.0143, -0.0039] (249/328 px, 9/15 blk) |
| Clay Loam | MAM | common | -0.0068 [-0.0168, +0.0006] (249/328 px, 9/15 blk) | -0.0017 [-0.0168, +0.0102] (249/328 px, 9/15 blk) | -0.0034 [-0.0117, +0.0028] (249/328 px, 9/15 blk) |
| Clay Loam | JJA | specific | +0.0109 [-0.0087, +0.0183] (96/156 px, 5/10 blk) | +0.0046 [+0.0000, +0.0100] (252/331 px, 9/15 blk) | +0.0010 [-0.0005, +0.0030] (249/328 px, 9/15 blk) |
| Clay Loam | JJA | common | +0.0108 [-0.0081, +0.0187] (95/153 px, 5/10 blk) | +0.0177 [-0.0082, +0.0276] (95/153 px, 5/10 blk) | -0.0052 [-0.0131, -0.0002] (95/153 px, 5/10 blk) |
| Clay Loam | SON | specific | -0.0024 [-0.0183, +0.0159] (252/325 px, 9/15 blk) | -0.0034 [-0.0071, +0.0027] (252/331 px, 9/15 blk) | -0.0022 [-0.0103, +0.0038] (249/328 px, 9/15 blk) |
| Clay Loam | SON | common | -0.0037 [-0.0179, +0.0138] (249/322 px, 9/15 blk) | -0.0149 [-0.0260, -0.0012] (249/322 px, 9/15 blk) | -0.0113 [-0.0219, +0.0021] (249/322 px, 9/15 blk) |
| Loam | All | specific | +0.0078 [-0.0112, +0.0228] (117/7467 px, 11/46 blk) | -0.0080 [-0.0142, -0.0023] (117/7467 px, 11/46 blk) | -0.0056 [-0.0112, -0.0009] (116/7265 px, 11/46 blk) |
| Loam | All | common | -0.0011 [-0.0112, +0.0063] (116/7235 px, 11/46 blk) | +0.0010 [-0.0048, +0.0058] (116/7235 px, 11/46 blk) | +0.0020 [-0.0029, +0.0060] (116/7235 px, 11/46 blk) |
| Loam | DJF | specific | +0.0104 [-0.0032, +0.0207] (117/6627 px, 11/44 blk) | -0.0169 [-0.0278, -0.0070] (117/7463 px, 11/46 blk) | -0.0082 [-0.0147, -0.0023] (116/7261 px, 11/46 blk) |
| Loam | DJF | common | +0.0108 [-0.0016, +0.0216] (113/6125 px, 11/44 blk) | +0.0011 [-0.0053, +0.0067] (113/6125 px, 11/44 blk) | +0.0011 [-0.0042, +0.0056] (113/6125 px, 11/44 blk) |
| Loam | MAM | specific | +0.0023 [-0.0228, +0.0196] (117/7447 px, 11/46 blk) | -0.0064 [-0.0137, +0.0016] (117/7467 px, 11/46 blk) | -0.0090 [-0.0161, -0.0023] (113/6916 px, 11/46 blk) |
| Loam | MAM | common | -0.0058 [-0.0218, +0.0088] (113/6693 px, 11/46 blk) | -0.0035 [-0.0162, +0.0121] (113/6693 px, 11/46 blk) | -0.0049 [-0.0116, +0.0046] (113/6693 px, 11/46 blk) |
| Loam | JJA | specific | +0.0145 [-0.0005, +0.0216] (95/6773 px, 11/46 blk) | +0.0027 [-0.0060, +0.0105] (99/7390 px, 11/46 blk) | -0.0012 [-0.0102, +0.0073] (97/7227 px, 11/46 blk) |
| Loam | JJA | common | -0.0501 [-0.0933, +0.0172] (56/6361 px, 11/44 blk) | +0.0074 [-0.0050, +0.0264] (56/6361 px, 11/44 blk) | +0.0008 [-0.0114, +0.0184] (56/6361 px, 11/44 blk) |
| Loam | SON | specific | +0.0252 [+0.0081, +0.0389] (117/7442 px, 11/46 blk) | +0.0019 [-0.0087, +0.0100] (117/7467 px, 11/46 blk) | +0.0111 [+0.0011, +0.0215] (116/7219 px, 11/46 blk) |
| Loam | SON | common | +0.0135 [-0.0011, +0.0239] (116/6871 px, 11/46 blk) | +0.0071 [-0.0000, +0.0123] (116/6871 px, 11/46 blk) | +0.0160 [+0.0040, +0.0232] (116/6871 px, 11/46 blk) |

---

## 11. Diagnostic du Dessèchement Estival (Figure F5)

Part des valeurs sous le seuil $\theta_{\text{wlt}} + 0.02$ calculée par mois pour le modèle Noah-MP Open Loop (`ssm_ol_daily`) et pour SMAP brut (`smap_sm`) sur la série 2016–2020 :

| Mois | Noah-MP Open Loop (% sous seuil) | SMAP brut (% sous seuil) |
|:---:|:---:|:---:|
| 01 | 6.46% | 36.86% |
| 02 | 8.98% | 43.26% |
| 03 | 12.75% | 47.27% |
| 04 | 13.06% | 49.73% |
| 05 | 19.74% | 56.56% |
| 06 | 43.19% | 62.98% |
| 07 | 79.80% | 58.16% |
| 08 | 62.52% | 50.05% |
| 09 | 26.62% | 47.77% |
| 10 | 18.93% | 45.30% |
| 11 | 12.74% | 39.36% |
| 12 | 6.84% | 26.97% |

- Figure associée : [F5_wilt_diagnostic_ol_vs_smap.png](file:///home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures/F5_wilt_diagnostic_ol_vs_smap.png) (version anglaise : [F5_wilt_diagnostic_ol_vs_smap_en.png](file:///home/mohammad.elaabaribao/.gemini/antigravity-ide/brain/a99d4e7d-be16-47d8-8e23-261668b8ac00/figures/en/F5_wilt_diagnostic_ol_vs_smap_en.png)).
- Tableau de données mensuel : [figures/descriptive/F5_wilt_diagnostic_monthly.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F5_wilt_diagnostic_monthly.csv).
- Données spatiales JJA : [figures/descriptive/F5_wilt_diagnostic_jja_map.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F5_wilt_diagnostic_jja_map.csv).

---

## 12. Effectifs de la Figure F3 (Nombre d'Observations et Pixels Distincts)

Extrait des effectifs sous la figure F3 issus de [figures/descriptive/F3_d_by_gmia_texture_season.csv](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/descriptive/F3_d_by_gmia_texture_season.csv) :

| Texture | Période | Classe GMIA | Nombre d'observations | Pixels distincts | Moyenne ± demi-IC 95% |
|:---|:---|:---:|:---:|:---:|:---:|
| Clay | Autres saisons (DJF, MAM, SON) | 0% | 146,401 | 136 | -0.0746 ± 0.0252 |
| Clay | Autres saisons (DJF, MAM, SON) | 0-10% | 266,655 | 240 | -0.0771 ± 0.0148 |
| Clay | Autres saisons (DJF, MAM, SON) | 10-30% | 64,222 | 55 | -0.0907 ± 0.0158 |
| Clay | Autres saisons (DJF, MAM, SON) | >30% | 29,506 | 26 | -0.1144 ± 0.0147 |
| Clay Loam | Autres saisons (DJF, MAM, SON) | 0% | 356,638 | 332 | -0.0739 ± 0.0090 |
| Clay Loam | Autres saisons (DJF, MAM, SON) | 0-10% | 400,541 | 356 | -0.0758 ± 0.0068 |
| Clay Loam | Autres saisons (DJF, MAM, SON) | 10-30% | 142,486 | 121 | -0.0794 ± 0.0063 |
| Clay Loam | Autres saisons (DJF, MAM, SON) | >30% | 309,922 | 252 | -0.0675 ± 0.0086 |
| Loam | Autres saisons (DJF, MAM, SON) | 0% | 8,008,936 | 7,480 | -0.0517 ± 0.0058 |
| Loam | Autres saisons (DJF, MAM, SON) | 0-10% | 3,087,116 | 3,151 | -0.0461 ± 0.0096 |
| Loam | Autres saisons (DJF, MAM, SON) | 10-30% | 330,828 | 319 | -0.0508 ± 0.0097 |
| Loam | Autres saisons (DJF, MAM, SON) | >30% | 121,308 | 117 | -0.0483 ± 0.0130 |
| Sandy Clay Loam | Autres saisons (DJF, MAM, SON) | 0% | 660,498 | 617 | -0.0104 ± 0.0094 |
| Sandy Clay Loam | Autres saisons (DJF, MAM, SON) | 0-10% | 68,025 | 66 | -0.0151 ± 0.0144 |
| Sandy Loam | Autres saisons (DJF, MAM, SON) | 0% | 351,867 | 312 | -0.0706 ± 0.0102 |
| Sandy Loam | Autres saisons (DJF, MAM, SON) | 0-10% | 80,957 | 74 | -0.0560 ± 0.0102 |
| Sandy Loam | Autres saisons (DJF, MAM, SON) | 10-30% | 11,297 | 10 | -0.0399 ± 0.0075 |
| Sandy Loam | Autres saisons (DJF, MAM, SON) | >30% | 16,974 | 14 | -0.0264 ± 0.0173 |
| Clay | JJA (Été) | 0% | 44,645 | 136 | -0.0971 ± 0.0227 |
| Clay | JJA (Été) | 0-10% | 86,156 | 240 | -0.0988 ± 0.0152 |
| Clay | JJA (Été) | 10-30% | 21,231 | 55 | -0.1047 ± 0.0128 |
| Clay | JJA (Été) | >30% | 9,274 | 26 | -0.1209 ± 0.0089 |
| Clay Loam | JJA (Été) | 0% | 112,331 | 332 | -0.0555 ± 0.0088 |
| Clay Loam | JJA (Été) | 0-10% | 127,043 | 356 | -0.0564 ± 0.0099 |
| Clay Loam | JJA (Été) | 10-30% | 44,404 | 121 | -0.0570 ± 0.0081 |
