# Rapport d'Évaluation Phase 2b : Premier Modèle d'Apprentissage (Random Forest)

**Date d'exécution** : 2026-11-23 12:54:29  
**Période d'entraînement** : 2016–2018 (`split == 'train'`).  
**Période d'évaluation** : 2019 (`split == 'val'`).  
**Statut de l'année 2020** : 2020 reste strictement scellée.  
**Sous-échantillonnage d'entraînement** : Échantillonnage stratifié uniforme de 20 observations par pixel sur les observations conformes aux règles QC hors flétrissement. Effectif d'entraînement = 285,000 observations réparties uniformément sur 14,250 mailles continentales.  
**Valeur de base (phi_0)** : -0.043058 m³/m³.

---

## 1. Validation Croisée Spatiale (5 Plis Équilibrés de 1°)

| Pli de validation | Effectif train | Effectif val | RMSE hors-pli | MAE hors-pli | Biais hors-pli | R² hors-pli |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 226,460 | 58,540 | 0.0245 | 0.0166 | -0.0019 | 0.8149 |
| 1 | 228,400 | 56,600 | 0.0279 | 0.0187 | -0.0050 | 0.7775 |
| 2 | 230,160 | 54,840 | 0.0271 | 0.0190 | -0.0004 | 0.7434 |
| 3 | 227,080 | 57,920 | 0.0256 | 0.0174 | +0.0020 | 0.7751 |
| 4 | 227,900 | 57,100 | 0.0256 | 0.0181 | +0.0056 | 0.7696 |

---

## 2. Attributions SHAP par Groupe et Stabilité Spatiale (Jalon H4)

Attributions SHAP moyennes calculées sur les échantillons de validation de chaque pli :

| Pli spatial | Groupe Removable (à retirer) | Groupe Retained (à conserver) | Groupe Context (analyse sensibilité) |
|:---|:---:|:---:|:---:|
| Plis 0 | 0.0139 (26.4%) | 0.0028 (5.4%) | 0.0361 (68.3%) |
| Plis 1 | 0.0141 (27.9%) | 0.0037 (7.3%) | 0.0329 (64.8%) |
| Plis 2 | 0.0127 (27.2%) | 0.0032 (6.8%) | 0.0307 (66.0%) |
| Plis 3 | 0.0121 (24.4%) | 0.0035 (6.9%) | 0.0342 (68.7%) |
| Plis 4 | 0.0091 (19.8%) | 0.0030 (6.5%) | 0.0338 (73.7%) |
| **Moyenne 5 plis** | **25.1%** | **6.6%** | **68.3%** |

---

## 3. Population d'Évaluation 2019 et Observations Assimilables Après Recalcul du Flétrissement

Socle d'observations 2019 passant les règles QC hors flétrissement : 3,660,799 observations (95.95% des 3,815,354 lignes brutes).

| Échantillon | Méthode | Observations assimilables | Part du socle hors flétrissement (3,660,799 obs) | Part du total 2019 (3,815,354 obs) |
|:---|:---|:---:|:---:|:---:|
| Spécifique | **B0 (NoCDF)** | 1,912,072 | 52.23% | 50.12% |
| Spécifique | **B1 (Climatologie)** | 3,133,887 | 85.61% | 82.14% |
| Spécifique | **B2 (CDF v2)** | 2,622,380 | 71.63% | 68.73% |
| Spécifique | **RF (SAR main)** | 3,548,975 | 96.95% | 93.02% |
| Spécifique | **RF-ctx (Sensibilité)** | 2,902,278 | 79.28% | 76.07% |
| Commun | **Commun (5 méthodes)** | 1,584,418 | 43.28% | 41.53% |

---

## 4. Métriques Globales Comparatives : Biais et RMSE Résiduels

| Méthode | Biais résiduel (Échantillon spécifique) | RMSE résiduelle (Échantillon spécifique) | Biais résiduel (Échantillon commun) | RMSE résiduelle (Échantillon commun) |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | -0.0265 m³/m³ | 0.0597 m³/m³ | -0.0341 m³/m³ | 0.0627 m³/m³ |
| **B1 (Climatologie)** | +0.0135 m³/m³ | 0.0427 m³/m³ | +0.0122 m³/m³ | 0.0441 m³/m³ |
| **B2 (CDF v2)** | +0.0086 m³/m³ | 0.0354 m³/m³ | +0.0096 m³/m³ | 0.0335 m³/m³ |
| **RF (SAR main)** | +0.0061 m³/m³ | 0.0456 m³/m³ | +0.0028 m³/m³ | 0.0478 m³/m³ |
| **RF-ctx (Sensibilité)** | +0.0062 m³/m³ | 0.0309 m³/m³ | +0.0147 m³/m³ | 0.0371 m³/m³ |

---

## 5. Performance Saisonnière sur l'Échantillon Commun (RMSE Résiduelle)

| Saison | RMSE B0 (NoCDF) | RMSE B1 (Climatologie) | RMSE B2 (CDF v2) | RMSE RF (SAR main) | RMSE RF-ctx (Sensibilité) | Observations communes |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| DJF | 0.0693 | 0.0408 | 0.0311 | 0.0460 | 0.0369 | 567,862 |
| MAM | 0.0662 | 0.0496 | 0.0318 | 0.0452 | 0.0380 | 430,612 |
| JJA | 0.0398 | 0.0335 | 0.0347 | 0.0555 | 0.0265 | 117,960 |
| SON | 0.0554 | 0.0450 | 0.0373 | 0.0503 | 0.0388 | 467,984 |

---

## 6. Corrélation d'Anomalies avec SMAP Brut ($r_{\text{anom}}$)

| Méthode | $r_{\text{anom}}$ spécifique | N spécifique | $r_{\text{anom}}$ commun | N commun |
|:---|:---:|:---:|:---:|:---:|
| **B0 (NoCDF)** | 1.0000 | 1,912,072 | 1.0000 | 1,584,418 |
| **B1 (Climatologie)** | 1.0000 | 3,133,887 | 1.0000 | 1,584,418 |
| **B2 (CDF v2)** | 0.8009 | 2,622,380 | 0.8195 | 1,584,418 |
| **RF (SAR main)** | 0.9869 | 3,548,975 | 0.9879 | 1,584,418 |
| **RF-ctx (Sensibilité)** | 0.9171 | 2,902,278 | 0.9354 | 1,584,418 |

---

## 7. Préservation du Signal d'Irrigation : Contraste GMIA >30 % vs 0 % en Clay Loam

Contraste $\Delta d_{\text{res}} = \bar{d}_{\text{res}}(>30\%) - \bar{d}_{\text{res}}(0\%)$ avec intervalle de confiance à 95% par bootstrap par blocs de 1° (500 itérations) :

| Saison | Échantillon | Contraste B0 | Contraste B1 | Contraste B2 | Contraste RF (SAR) | Effectif (>30% / 0%) |
|:---|:---|:---:|:---:|:---:|:---:|:---:|
| All | specific | -0.0025 [-0.0169, +0.0103] | -0.0049 [-0.0076, -0.0019] | -0.0063 [-0.0123, -0.0027] | +0.0054 [-0.0037, +0.0109] | 252/331 px, 9/15 blk |
| All | common | -0.0034 [-0.0175, +0.0079] | -0.0128 [-0.0212, -0.0042] | -0.0094 [-0.0195, +0.0014] | +0.0029 [-0.0104, +0.0139] | 249/328 px, 9/15 blk |
| DJF | specific | -0.0010 [-0.0219, +0.0167] | -0.0105 [-0.0206, -0.0028] | -0.0074 [-0.0161, -0.0005] | +0.0080 [-0.0081, +0.0210] | 252/331 px, 9/15 blk |
| DJF | common | -0.0020 [-0.0218, +0.0145] | -0.0164 [-0.0310, -0.0050] | -0.0105 [-0.0231, +0.0028] | +0.0031 [-0.0147, +0.0182] | 249/328 px, 9/15 blk |
| MAM | specific | -0.0062 [-0.0158, +0.0014] | -0.0032 [-0.0125, +0.0056] | -0.0096 [-0.0151, -0.0021] | +0.0016 [-0.0038, +0.0055] | 252/331 px, 9/15 blk |
| MAM | common | -0.0068 [-0.0154, +0.0006] | -0.0017 [-0.0170, +0.0098] | -0.0034 [-0.0125, +0.0031] | +0.0042 [-0.0056, +0.0098] | 249/328 px, 9/15 blk |
| JJA | specific | +0.0109 [-0.0095, +0.0176] | +0.0046 [+0.0003, +0.0098] | +0.0010 [-0.0005, +0.0036] | +0.0036 [-0.0011, +0.0107] | 252/309 px, 9/15 blk |
| JJA | common | +0.0108 [-0.0097, +0.0178] | +0.0177 [-0.0091, +0.0271] | -0.0052 [-0.0138, -0.0004] | +0.0174 [-0.0090, +0.0245] | 95/153 px, 5/10 blk |
| SON | specific | -0.0024 [-0.0177, +0.0188] | -0.0034 [-0.0080, +0.0025] | -0.0022 [-0.0099, +0.0040] | +0.0039 [-0.0047, +0.0141] | 252/331 px, 9/15 blk |
| SON | common | -0.0037 [-0.0179, +0.0163] | -0.0149 [-0.0258, -0.0023] | -0.0113 [-0.0213, +0.0018] | +0.0020 [-0.0130, +0.0172] | 249/322 px, 9/15 blk |

---

## 8. Part des Valeurs Sous le Seuil de Flétrissement par Saison ($\text{SM}^* < \theta_{\text{wlt}} + 0.02$)

| Saison | B0 (NoCDF) | B1 (Climatologie) | B2 (CDF v2) | RF (SAR main) | RF-ctx (Sensibilité) |
|:---|:---:|:---:|:---:|:---:|
| DJF | 39.65% | 0.31% | 3.87% | 0.97% | 6.27% |
| MAM | 52.94% | 11.63% | 8.06% | 4.28% | 11.41% |
| JJA | 57.25% | 54.02% | 57.86% | 5.76% | 56.85% |
| SON | 44.46% | 5.91% | 16.17% | 2.14% | 21.60% |
## 9. Fichiers et Tables Produites

- **Table de métriques détaillées** : `reports/phase2b_rf_metrics.csv`
- **Table de stabilité SHAP 5 plis** : `reports/phase2b_rf_shap_stability.csv`
- **Table de contrastes GMIA** : `reports/phase2b_rf_gmia_contrast.csv`
