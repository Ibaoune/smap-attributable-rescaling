# Rapport de Phase 3 : Évaluation Scellée sur l'Année 2020

**Date d'exécution :** 2026-11-23 19:26:17  
**Période d'évaluation :** 2020-01-01 au 2020-12-31 (année 2020 scellée au préalable, ouverte une seule fois)  
**Modèle utilisé :** Modèle LightGBM pré-entraîné sur 2016–2018 (300 arbres, learning_rate=0,05, num_leaves=63, aucun réentraînement, aucun changement d'hyperparamètres)  
**Référence principale :** CDF B2-train reconstruite avec LDT strictement sur 2016-01-01 → 2019-01-01  
**Méthode candidate formelle (Décision 15) :** Variante V2 (Removable + Context Representation) : $\mathrm{SMAP}^*_{\mathrm{V2}} = \mathrm{SMAP} - (\Phi_{\mathrm{removable}} + \Phi_{\mathrm{context\_rep}} + \phi_0)$  
**Observations 2020 brutes :** 4,135,833  
**Observations 2020 passant QC (hors flétrissement) :** 3,976,522 (96.15%)  
**Constante de base :** $\phi_0 = -0.0430$ m³/m³ (TreeSHAP exact)

---

## 1. Évaluation des Critères de Succès Figés (Décision 15)

Verdict global avant ouverture : **CRITÈRE DE SUCCÈS DÉCISION 15 SATISFAIT**

| Critère | Formulation | Valeurs mesurées (V2 vs B2-train) | Résultat |
| :--- | :--- | :--- | :---: |
| 1. Biais absolu (échantillon spécifique) | $|\text{Biais}(\mathrm{V2})| < |\text{Biais}(\mathrm{B2\text{-}train})|$ | $|+0.0092| < |+0.0234|$ | **Validé** |
| 2. Biais absolu (échantillon commun) | $|\text{Biais}(\mathrm{V2})| < |\text{Biais}(\mathrm{B2\text{-}train})|$ | $|+0.0177| < |+0.0184|$ | **Validé** |
| 3. Corrélation d'anomalies (échantillon spécifique) | $r_{\mathrm{anom}}(\mathrm{V2}) > r_{\mathrm{anom}}(\mathrm{B2\text{-}train})$ | 0.8841 > 0.7893 | **Validé** |
| 4. Corrélation d'anomalies (échantillon commun) | $r_{\mathrm{anom}}(\mathrm{V2}) > r_{\mathrm{anom}}(\mathrm{B2\text{-}train})$ | 0.8981 > 0.8135 | **Validé** |
| 5. RMSE (échantillon commun) | $\text{RMSE}(\mathrm{V2}) \le 1{,}05 \times \text{RMSE}(\mathrm{B2\text{-}train})$ | 0.0375 m³/m³ $\le$ 0.0428 m³/m³ (B2-train: 0.0408) | **Validé** |

---

## 2. Métriques Complètes du Protocole 2a (Année 2020)

Le flétrissement ($\mathrm{SMAP}^* \ge \theta_{\mathrm{wlt}} + 0{,}02$) est recalculé sur l'observation corrigée pour chaque référence.

### 2.1 Échantillon assimilable spécifique à chaque référence

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 2,094,747 (52.7%) | 1,881,775 (47.3%) | -0.0179 m³/m³ | 0.0558 m³/m³ | 0.0528 m³/m³ | 0.6278 | 1.0000 | +0.0200 [-0.0011, +0.0263] (interprétable) |
| B1_Climatology | 3,226,787 (81.1%) | 749,735 (18.9%) | +0.0183 m³/m³ | 0.0454 m³/m³ | 0.0415 m³/m³ | 0.7546 | 1.0000 | +0.0136 [+0.0042, +0.0213] (interprétable) |
| B2_CDF_Legacy | 2,782,841 (70.0%) | 797,158 (20.0%) | +0.0161 m³/m³ | 0.0400 m³/m³ | 0.0367 m³/m³ | 0.7875 | 0.8046 | +0.0093 [+0.0016, +0.0141] (interprétable) |
| B2_train_CDF | 2,927,357 (73.6%) | 763,665 (19.2%) | +0.0234 m³/m³ | 0.0476 m³/m³ | 0.0415 m³/m³ | 0.7353 | 0.7893 | +0.0112 [+0.0023, +0.0168] (interprétable) |
| LGB_V1_Removable | 3,835,152 (96.4%) | 141,370 (3.6%) | +0.0178 m³/m³ | 0.0512 m³/m³ | 0.0480 m³/m³ | 0.6314 | 0.9677 | +0.0128 [+0.0033, +0.0205] (interprétable) |
| LGB_V2_Rem_ContextRep | 2,808,751 (70.6%) | 1,167,771 (29.4%) | +0.0092 m³/m³ | 0.0308 m³/m³ | 0.0294 m³/m³ | 0.8777 | 0.8841 | +0.0124 [+0.0064, +0.0182] (interprétable) |
| LGB_V3_Rem_AllContext | 2,719,931 (68.4%) | 1,256,591 (31.6%) | +0.0071 m³/m³ | 0.0287 m³/m³ | 0.0278 m³/m³ | 0.8953 | 0.8689 | +0.0061 [+0.0022, +0.0103] (interprétable) |

### 2.2 Échantillon commun (1,555,724 observations assimilables par toutes les références)

| Référence | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) | Contraste GMIA Clay Loam JJA |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | 1,555,724 (39.1%) | 0 (0.0%) | -0.0309 m³/m³ | 0.0601 m³/m³ | 0.0515 m³/m³ | 0.5647 | 1.0000 | +0.0197 [-0.0017, +0.0263] (interprétable) |
| B1_Climatology | 1,555,724 (39.1%) | 0 (0.0%) | +0.0156 m³/m³ | 0.0458 m³/m³ | 0.0431 m³/m³ | 0.7396 | 1.0000 | +0.0265 [+0.0039, +0.0372] (interprétable) |
| B2_CDF_Legacy | 1,555,724 (39.1%) | 0 (0.0%) | +0.0142 m³/m³ | 0.0356 m³/m³ | 0.0327 m³/m³ | 0.8067 | 0.8305 | +0.0121 [-0.0005, +0.0207] (interprétable) |
| B2_train_CDF | 1,555,724 (39.1%) | 0 (0.0%) | +0.0184 m³/m³ | 0.0408 m³/m³ | 0.0364 m³/m³ | 0.7650 | 0.8135 | +0.0147 [-0.0002, +0.0241] (interprétable) |
| LGB_V1_Removable | 1,555,724 (39.1%) | 0 (0.0%) | +0.0006 m³/m³ | 0.0471 m³/m³ | 0.0471 m³/m³ | 0.5785 | 0.9729 | +0.0282 [-0.0029, +0.0379] (interprétable) |
| LGB_V2_Rem_ContextRep | 1,555,724 (39.1%) | 0 (0.0%) | +0.0177 m³/m³ | 0.0375 m³/m³ | 0.0331 m³/m³ | 0.8517 | 0.8981 | +0.0129 [+0.0039, +0.0156] (interprétable) |
| LGB_V3_Rem_AllContext | 1,555,724 (39.1%) | 0 (0.0%) | +0.0140 m³/m³ | 0.0345 m³/m³ | 0.0315 m³/m³ | 0.8796 | 0.8794 | +0.0035 [-0.0053, +0.0092] (interprétable) |

---

## 3. Corrections Saisonnières Moyennes Appliquées $\overline{\Delta \mathrm{SMAP}}$ et Biais Résiduel $\overline{d_{\mathrm{res}}}$ (2020)

| Saison | Nombre d'obs | V1 : Correction | V1 : Biais $d_{\mathrm{res}}$ | V2 : Correction | V2 : Biais $d_{\mathrm{res}}$ | V3 : Correction | V3 : Biais $d_{\mathrm{res}}$ | B2-train : Correction | B2-train : Biais $d_{\mathrm{res}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| DJF | 977,929 | +0.0412 m³/m³ | -0.0023 m³/m³ | +0.0523 m³/m³ | +0.0089 m³/m³ | +0.0485 m³/m³ | +0.0051 m³/m³ | +0.0741 m³/m³ | +0.0284 m³/m³ |
| MAM | 989,504 | +0.0454 m³/m³ | -0.0028 m³/m³ | +0.0556 m³/m³ | +0.0074 m³/m³ | +0.0554 m³/m³ | +0.0072 m³/m³ | +0.0658 m³/m³ | +0.0156 m³/m³ |
| JJA | 1,045,144 | +0.0447 m³/m³ | +0.0405 m³/m³ | +0.0055 m³/m³ | +0.0014 m³/m³ | +0.0033 m³/m³ | -0.0009 m³/m³ | +0.0051 m³/m³ | +0.0004 m³/m³ |
| SON | 963,945 | +0.0471 m³/m³ | +0.0247 m³/m³ | +0.0305 m³/m³ | +0.0080 m³/m³ | +0.0270 m³/m³ | +0.0045 m³/m³ | +0.0474 m³/m³ | +0.0249 m³/m³ |
| Annual | 3,976,522 | +0.0446 m³/m³ | +0.0154 m³/m³ | +0.0356 m³/m³ | +0.0063 m³/m³ | +0.0331 m³/m³ | +0.0039 m³/m³ | +0.0471 m³/m³ | +0.0168 m³/m³ |

---

## 4. Stratification Saisonnière des Métriques (Échantillon spécifique)

| Référence | Saison | Observations assimilables | Rejets flétrissement | Biais ($d_{\mathrm{res}} = \mathrm{SMAP}^* - \mathrm{SSM}_{\mathrm{OL}}$) | RMSE | ubRMSE | $r_{\mathrm{SMAP}^* \mathrm{vs} \mathrm{OL}}$ | $r_{\mathrm{anom}}$ (SMAP brut) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| B0_Raw | DJF | 599,542 (61.3%) | 378,387 (38.7%) | -0.0353 m³/m³ | 0.0621 m³/m³ | 0.0511 m³/m³ | 0.6152 | 1.0000 |
| B0_Raw | MAM | 515,415 (52.1%) | 474,089 (47.9%) | -0.0373 m³/m³ | 0.0651 m³/m³ | 0.0533 m³/m³ | 0.6812 | 1.0000 |
| B0_Raw | JJA | 461,668 (44.2%) | 583,476 (55.8%) | +0.0169 m³/m³ | 0.0411 m³/m³ | 0.0374 m³/m³ | 0.5045 | 1.0000 |
| B0_Raw | SON | 518,122 (53.8%) | 445,823 (46.2%) | -0.0095 m³/m³ | 0.0490 m³/m³ | 0.0481 m³/m³ | 0.5849 | 1.0000 |
| B0_Raw | Annual | 2,094,747 (52.7%) | 1,881,775 (47.3%) | -0.0179 m³/m³ | 0.0558 m³/m³ | 0.0528 m³/m³ | 0.6278 | 1.0000 |
| B1_Climatology | DJF | 975,096 (99.7%) | 2,833 (0.3%) | +0.0247 m³/m³ | 0.0465 m³/m³ | 0.0393 m³/m³ | 0.7464 | 1.0000 |
| B1_Climatology | MAM | 886,230 (89.6%) | 103,274 (10.4%) | +0.0148 m³/m³ | 0.0483 m³/m³ | 0.0460 m³/m³ | 0.7607 | 1.0000 |
| B1_Climatology | JJA | 462,250 (44.2%) | 582,894 (55.8%) | +0.0029 m³/m³ | 0.0286 m³/m³ | 0.0285 m³/m³ | 0.7975 | 1.0000 |
| B1_Climatology | SON | 903,211 (93.7%) | 60,734 (6.3%) | +0.0227 m³/m³ | 0.0481 m³/m³ | 0.0424 m³/m³ | 0.6501 | 1.0000 |
| B1_Climatology | Annual | 3,226,787 (81.1%) | 749,735 (18.9%) | +0.0183 m³/m³ | 0.0454 m³/m³ | 0.0415 m³/m³ | 0.7546 | 1.0000 |
| B2_CDF_Legacy | DJF | 853,536 (87.3%) | 22,434 (2.3%) | +0.0205 m³/m³ | 0.0397 m³/m³ | 0.0341 m³/m³ | 0.7732 | 0.7802 |
| B2_CDF_Legacy | MAM | 787,847 (79.6%) | 48,982 (5.0%) | +0.0124 m³/m³ | 0.0407 m³/m³ | 0.0388 m³/m³ | 0.7827 | 0.8294 |
| B2_CDF_Legacy | JJA | 401,855 (38.4%) | 593,592 (56.8%) | +0.0053 m³/m³ | 0.0309 m³/m³ | 0.0304 m³/m³ | 0.7633 | 0.7756 |
| B2_CDF_Legacy | SON | 739,603 (76.7%) | 132,150 (13.7%) | +0.0208 m³/m³ | 0.0439 m³/m³ | 0.0386 m³/m³ | 0.7129 | 0.8070 |
| B2_CDF_Legacy | Annual | 2,782,841 (70.0%) | 797,158 (20.0%) | +0.0161 m³/m³ | 0.0400 m³/m³ | 0.0367 m³/m³ | 0.7875 | 0.8046 |
| B2_train_CDF | DJF | 883,838 (90.4%) | 12,811 (1.3%) | +0.0292 m³/m³ | 0.0482 m³/m³ | 0.0383 m³/m³ | 0.7203 | 0.7836 |
| B2_train_CDF | MAM | 851,133 (86.0%) | 70,110 (7.1%) | +0.0191 m³/m³ | 0.0484 m³/m³ | 0.0444 m³/m³ | 0.7365 | 0.8093 |
| B2_train_CDF | JJA | 416,752 (39.9%) | 579,133 (55.4%) | +0.0063 m³/m³ | 0.0332 m³/m³ | 0.0326 m³/m³ | 0.7309 | 0.7439 |
| B2_train_CDF | SON | 775,634 (80.5%) | 101,611 (10.5%) | +0.0307 m³/m³ | 0.0525 m³/m³ | 0.0426 m³/m³ | 0.6437 | 0.7958 |
| B2_train_CDF | Annual | 2,927,357 (73.6%) | 763,665 (19.2%) | +0.0234 m³/m³ | 0.0476 m³/m³ | 0.0415 m³/m³ | 0.7353 | 0.7893 |
| LGB_V1_Removable | DJF | 962,377 (98.4%) | 15,552 (1.6%) | -0.0012 m³/m³ | 0.0445 m³/m³ | 0.0445 m³/m³ | 0.6528 | 0.9579 |
| LGB_V1_Removable | MAM | 953,407 (96.4%) | 36,097 (3.6%) | -0.0003 m³/m³ | 0.0515 m³/m³ | 0.0515 m³/m³ | 0.6485 | 0.9779 |
| LGB_V1_Removable | JJA | 984,840 (94.2%) | 60,304 (5.8%) | +0.0454 m³/m³ | 0.0577 m³/m³ | 0.0356 m³/m³ | 0.3854 | 0.9080 |
| LGB_V1_Removable | SON | 934,528 (96.9%) | 29,417 (3.1%) | +0.0269 m³/m³ | 0.0501 m³/m³ | 0.0423 m³/m³ | 0.6082 | 0.9713 |
| LGB_V1_Removable | Annual | 3,835,152 (96.4%) | 141,370 (3.6%) | +0.0178 m³/m³ | 0.0512 m³/m³ | 0.0480 m³/m³ | 0.6314 | 0.9677 |
| LGB_V2_Rem_ContextRep | DJF | 842,027 (86.1%) | 135,902 (13.9%) | +0.0102 m³/m³ | 0.0322 m³/m³ | 0.0306 m³/m³ | 0.8404 | 0.9033 |
| LGB_V2_Rem_ContextRep | MAM | 814,327 (82.3%) | 175,177 (17.7%) | +0.0087 m³/m³ | 0.0342 m³/m³ | 0.0331 m³/m³ | 0.8694 | 0.9017 |
| LGB_V2_Rem_ContextRep | JJA | 440,551 (42.2%) | 604,593 (57.8%) | +0.0049 m³/m³ | 0.0201 m³/m³ | 0.0195 m³/m³ | 0.8975 | 0.6452 |
| LGB_V2_Rem_ContextRep | SON | 711,846 (73.8%) | 252,099 (26.2%) | +0.0113 m³/m³ | 0.0304 m³/m³ | 0.0282 m³/m³ | 0.8690 | 0.8823 |
| LGB_V2_Rem_ContextRep | Annual | 2,808,751 (70.6%) | 1,167,771 (29.4%) | +0.0092 m³/m³ | 0.0308 m³/m³ | 0.0294 m³/m³ | 0.8777 | 0.8841 |
| LGB_V3_Rem_AllContext | DJF | 830,846 (85.0%) | 147,083 (15.0%) | +0.0064 m³/m³ | 0.0282 m³/m³ | 0.0275 m³/m³ | 0.8735 | 0.8808 |
| LGB_V3_Rem_AllContext | MAM | 800,363 (80.9%) | 189,141 (19.1%) | +0.0092 m³/m³ | 0.0339 m³/m³ | 0.0327 m³/m³ | 0.8867 | 0.8891 |
| LGB_V3_Rem_AllContext | JJA | 406,649 (38.9%) | 638,495 (61.1%) | +0.0033 m³/m³ | 0.0180 m³/m³ | 0.0177 m³/m³ | 0.9059 | 0.6213 |
| LGB_V3_Rem_AllContext | SON | 682,073 (70.8%) | 281,872 (29.2%) | +0.0080 m³/m³ | 0.0276 m³/m³ | 0.0265 m³/m³ | 0.8889 | 0.8692 |
| LGB_V3_Rem_AllContext | Annual | 2,719,931 (68.4%) | 1,256,591 (31.6%) | +0.0071 m³/m³ | 0.0287 m³/m³ | 0.0278 m³/m³ | 0.8953 | 0.8689 |

---

## 5. Figures de Publication (Année 2020)

- **Cartes saisonnières de correction moyenne (V1, V2, V3, B2-train) :** [F11_test2020_seasonal_applied_corrections.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/F11_test2020_seasonal_applied_corrections.png) | EN: [F11_test2020_seasonal_applied_corrections_en.png](file:///srv/lustre01/project/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/smap-attributable-rescaling/figures/evaluation/en/F11_test2020_seasonal_applied_corrections_en.png)
