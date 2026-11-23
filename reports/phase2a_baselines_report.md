# Rapport de Synthèse Phase 2a : Évaluation des Références (Baselines B0, B1, B2)

**Date d'exécution** : 2026-11-23 00:04:49  
**Protocole de scellement temporel** :
- **Apprentissage (train)** : 2016–2018 (`split == 'train'`).
- **Évaluation (val)** : Année 2019 (`split == 'val'`).
- **Test scellé (test)** : Année 2020 reste **strictement scellée** et intacte.
- **Règle absolue** : Aucun chiffre écrit à la main. Tous les indicateurs sont calculés dynamiquement à partir des données réelles.

---

## 1. Description des Références Évaluées

1. **B0 (NoCDF - Pas de correction)** :
   - On assimile $d$ brut tel quel : l'observation corrigée est $\text{smap}^* = \text{smap_sm}$, et le résidu est $d_{\text{res}} = d = \text{smap_sm} - \text{ssm_ol_daily}$.
2. **B1 (Climatologie pixel-mois)** :
   - Moyenne locale de $d$ par pixel et par mois, $\hat{c}(p, m)$, apprise exclusivement sur 2016–2018 (filtrée par le contrôle qualité LIS).
   - Sur 2019 : $\text{smap}^* = \text{smap_sm} - \hat{c}(p, m)$ et $d_{\text{res}} = d - \hat{c}(p, m)$.
3. **B2 (CDF matching v2 appliqué hors ligne)** :
   - Reproduction exacte bit-à-bit de la sous-routine Fortran `LIS_rescale_with_CDF_matching` avec les fichiers LDT officiels (`noahmp_model_cdf_NorthMor.nc` et `smap_obs_cdf_NorthMor.nc`).
   - Sur 2019 : $\text{smap}^* = \text{CDF_rescale}(\text{smap_sm}, p, m)$ et $d_{\text{res}} = \text{smap}^* - \text{ssm_ol_daily}$.

---

## 2. Vérification Bit-à-Bit du CDF Matching contre DAOBS LIS (Record 2)

La transformation CDF matching hors ligne a été rigoureusement confrontée à l'enregistrement 2 des fichiers binaires DAOBS issus de la simulation LIS opérationnelle `DA_cdf_noirr_2016_2020` sur un échantillon multi-saisonnier de 5 dates :

| Date et Heure DAOBS | Nombre de points comparés | Biais absolu moyen | Différence absolue max | Concordance exacte (< 1e-5) |
|:---|:---:|:---:|:---:|:---:|
| `2016-01-01 06:00` | 5,420 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-01-01 06:00` | 4,011 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-04-01 06:00` | 8,532 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-07-24 06:00` | 398 | 0.00000000 | 0.00000000 | 100.00% |
| `2019-10-01 18:00` | 647 | 0.00000000 | 0.00000000 | 100.00% |

> **Conclusion de validation B2** : L'implémentation hors ligne reproduit fidèlement et à l'identique les observations remises à l'échelle injectées dans le filtre d'assimilation EnKF de LIS.

---

## 3. Résultats Globaux sur l'Année d'Évaluation 2019

Sur l'ensemble des 1,912,072 observations assimilables (après filtres de qualité LIS) de l'année 2019 :

| Référence | Biais Résiduel ($d_{\text{res}}$) | RMSE Résiduelle ($d_{\text{res}}$) | Part sous $\theta_{\text{wlt}} + 0.02$ (Toutes obs) | Part sous $\theta_{\text{wlt}} + 0.02$ (Obs QC) |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | -0.0265 m³/m³ | 0.0597 m³/m³ | 46.40% | 0.00% |
| **B1 (Climatologie)** | +0.0079 m³/m³ | 0.0416 m³/m³ | 18.99% | 8.93% |
| **B2 (CDF v2 hors ligne)** | +0.0079 m³/m³ | 0.0329 m³/m³ | 17.37% | 8.76% |

### Analyse des Performances Globales
- **B0 (NoCDF)** conserve un biais négatif persistant (-0.0265 m³/m³) et une RMSE de 0.0597 m³/m³.
- **B1 (Climatologie pixel-mois)** élimine quasi-parfaitement le biais global moyen (+0.0079 m³/m³) et réduit la RMSE à 0.0416 m³/m³.
- **B2 (CDF v2 matching)** ramène le biais moyen à +0.0079 m³/m³ et la RMSE à 0.0329 m³/m³. Cependant, elle projette 17.37% de l'ensemble des observations sous le seuil critique de flétrissement $\theta_{\text{wlt}} + 0.02$, expliquant le rejet massif d'observations par LIS.

---

## 4. Résultats Stratifiés par Saison (Année 2019)

| Saison | Observations QC | RMSE B0 | RMSE B1 | RMSE B2 | Rejet Flétrissement B0 | Rejet Flétrissement B1 | Rejet Flétrissement B2 |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| DJF | 600,414 | 0.0684 | 0.0413 | 0.0314 | 38.30% | 3.06% | 3.74% |
| MAM | 484,559 | 0.0645 | 0.0465 | 0.0320 | 51.02% | 16.38% | 7.71% |
| JJA | 273,124 | 0.0388 | 0.0279 | 0.0285 | 56.00% | 60.68% | 56.60% |
| SON | 553,975 | 0.0532 | 0.0429 | 0.0370 | 43.48% | 10.73% | 15.77% |

---

## 5. Résultats Stratifiés par Texture Dominante (Top 5 Classes)

| Texture du Sol | Observations QC | Biais B0 | RMSE B0 | Biais B1 | RMSE B1 | Biais B2 | RMSE B2 | Part < $\theta_{\text{wlt}} + 0.02$ (B2) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Loam | 1,430,068 | -0.0291 | 0.0613 | +0.0083 | 0.0418 | +0.0083 | 0.0334 | 19.03% |
| Clay Loam | 98,285 | -0.0510 | 0.0702 | +0.0070 | 0.0437 | +0.0073 | 0.0350 | 10.65% |
| Sandy Clay Loam | 137,156 | +0.0028 | 0.0449 | +0.0099 | 0.0416 | +0.0067 | 0.0281 | 19.36% |
| Sandy Loam | 53,757 | -0.0440 | 0.0699 | +0.0146 | 0.0482 | +0.0146 | 0.0338 | 13.02% |
| Clay | 46,540 | -0.0433 | 0.0658 | +0.0104 | 0.0454 | +0.0010 | 0.0333 | 0.00% |

---

## 6. Résultats Stratifiés par Classe d'Irrigation GMIA

| Classe GMIA | Observations QC | Biais B0 | RMSE B0 | Biais B1 | RMSE B1 | Biais B2 | RMSE B2 | Part < $\theta_{\text{wlt}} + 0.02$ (B2) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0% | 1,211,241 | -0.0239 | 0.0578 | +0.0089 | 0.0419 | +0.0092 | 0.0333 | 18.44% |
| 0-10% | 601,277 | -0.0298 | 0.0622 | +0.0062 | 0.0408 | +0.0055 | 0.0315 | 13.86% |
| 10-30% | 56,523 | -0.0381 | 0.0667 | +0.0067 | 0.0426 | +0.0080 | 0.0356 | 19.73% |
| >30% | 43,031 | -0.0387 | 0.0652 | +0.0048 | 0.0420 | +0.0059 | 0.0369 | 20.82% |

---

## 7. Tableau Exhaustif de Toutes les Strate Évaluées

| Type de Strate | Nom de la Strate | N Obs QC | Biais B0 | RMSE B0 | Flétr. B0 | Biais B1 | RMSE B1 | Flétr. B1 | Biais B2 | RMSE B2 | Flétr. B2 |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Overall | Global (All) | 1,912,072 | -0.0265 | 0.0597 | 46.40% | +0.0079 | 0.0416 | 18.99% | +0.0079 | 0.0329 | 17.37% |
| Season | DJF | 600,414 | -0.0449 | 0.0684 | 38.30% | +0.0083 | 0.0413 | 3.06% | +0.0058 | 0.0314 | 3.74% |
| Season | MAM | 484,559 | -0.0441 | 0.0645 | 51.02% | +0.0029 | 0.0465 | 16.38% | +0.0083 | 0.0320 | 7.71% |
| Season | JJA | 273,124 | +0.0199 | 0.0388 | 56.00% | +0.0008 | 0.0279 | 60.68% | +0.0053 | 0.0285 | 56.60% |
| Season | SON | 553,975 | -0.0140 | 0.0532 | 43.48% | +0.0152 | 0.0429 | 10.73% | +0.0112 | 0.0370 | 15.77% |
| Texture | Loam | 1,430,068 | -0.0291 | 0.0613 | 47.15% | +0.0083 | 0.0418 | 18.57% | +0.0083 | 0.0334 | 19.03% |
| Texture | Clay Loam | 98,285 | -0.0510 | 0.0702 | 64.02% | +0.0070 | 0.0437 | 23.39% | +0.0073 | 0.0350 | 10.65% |
| Texture | Sandy Clay Loam | 137,156 | +0.0028 | 0.0449 | 21.56% | +0.0099 | 0.0416 | 16.44% | +0.0067 | 0.0281 | 19.36% |
| Texture | Sandy Loam | 53,757 | -0.0440 | 0.0699 | 50.71% | +0.0146 | 0.0482 | 12.35% | +0.0146 | 0.0338 | 13.02% |
| Texture | Clay | 46,540 | -0.0433 | 0.0658 | 61.37% | +0.0104 | 0.0454 | 24.60% | +0.0010 | 0.0333 | 0.00% |
| GMIA | 0% | 1,211,241 | -0.0239 | 0.0578 | 49.26% | +0.0089 | 0.0419 | 19.85% | +0.0092 | 0.0333 | 18.44% |
| GMIA | 0-10% | 601,277 | -0.0298 | 0.0622 | 35.47% | +0.0062 | 0.0408 | 15.01% | +0.0055 | 0.0315 | 13.86% |
| GMIA | 10-30% | 56,523 | -0.0381 | 0.0667 | 56.86% | +0.0067 | 0.0426 | 22.93% | +0.0080 | 0.0356 | 19.73% |
| GMIA | >30% | 43,031 | -0.0387 | 0.0652 | 63.56% | +0.0048 | 0.0420 | 28.77% | +0.0059 | 0.0369 | 20.82% |

---

## 8. Synthèse pour le Point de Contrôle 5

1. **Reproduction exacte de B2** : Le module CDF matching hors ligne est validé à 100% contre les enregistrements DAOBS de LIS.
2. **Benchmark établi** : Les baselines B0, B1 et B2 fournissent les seuils de référence indispensables pour jauger le gain apporté par les futurs modèles de ML explicables (LightGBM / TabNet).
3. **Maintien de l'étanchéité** : L'année 2020 demeure strictement scellée.
