# Rapport de Diagnostic Phase 2a-bis : Signal d'Irrigation Dynamique

**Période d'évaluation** : 2016-01-01 au 2019-12-31 (4 années complètes).  
**Statut de l'année 2020** : 2020 reste strictement scellée.  
**Population d'observation** : Paires d'observations SMAP successives du même passage orbital (AM avec AM, PM avec PM), espacées de 1 à 3 jours, issues des observations conformes aux règles QC hors flétrissement.  
**Effectif total analysé** : 13,878,569 paires successives (14,250 mailles distinctes, 50 blocs spatiaux distincts de 1°).

---

## 1. Définition Opérationnelle et Seuils

Un événement d'humidification sans précipitation est défini par :
- Un intervalle temporel de 1 à 3 jours entre deux observations du même passage au même pixel.
- Une élévation d'humidité du sol : $\Delta 	ext{SM} > 	au_{	ext{sm}}$.
- Un cumul de précipitations IMERG sur l'intervalle : $	ext{Rain}_{t_1 	o t_2} < 	au_{	ext{rain}}$.
- Configuration par défaut : $	au_{	ext{sm}} = 0.03	ext{ m}^3/	ext{m}^3$, $	au_{	ext{rain}} = 1.0	ext{ mm}$.
- Le modèle Open Loop (Noah-MP sans irrigation) subit exactement le même calcul sur les mêmes paires d'observations pour servir de témoin négatif.

---

## 2. Statistiques Globales (Configuration par Défaut)

- **Nombre total de paires de passages successifs** : 13,878,569
- **Événements détectés sur SMAP (ensemble du domaine et de l'année)** : 261,486 (1.884%, IC 95% : [1.512%, 2.326%])
- **Événements détectés sur l'Open Loop (témoin)** : 0 (0.000%, IC 95% : [0.000%, 0.000%])
- **Écart global SMAP − OL** : +1.884%

---

## 3. Signal d'Été (JJA) par Classe d'Irrigation GMIA à Texture Fixée

### 3.1. Texture Clay Loam (JJA)

| Classe GMIA | Paires | Pixels | Blocs 1° | Événements SMAP | Fréquence SMAP [IC 95%] | Événements OL | Fréquence OL [IC 95%] | Contraste SMAP − OL |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0% | 76,759 | 329 | 15 | 1,614 | 2.103% [1.187%, 3.213%] | 0 | 0.000% [0.000%, 0.000%] | +2.103% |
| 0-10% | 88,753 | 354 | 17 | 2,018 | 2.274% [1.083%, 3.949%] | 0 | 0.000% [0.000%, 0.000%] | +2.274% |
| 10-30% | 31,458 | 121 | 10 | 1,217 | 3.869% [0.881%, 9.843%] | 0 | 0.000% [0.000%, 0.000%] | +3.869% |
| >30% | 70,316 | 252 | 9 | 3,186 | 4.531% [1.015%, 12.572%] | 0 | 0.000% [0.000%, 0.000%] | +4.531% |

### 3.2. Texture Loam (JJA)

| Classe GMIA | Paires | Pixels | Blocs 1° | Événements SMAP | Fréquence SMAP [IC 95%] | Événements OL | Fréquence OL [IC 95%] | Contraste SMAP − OL |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0% | 1,796,821 | 7,467 | 46 | 19,915 | 1.108% [0.793%, 1.500%] | 0 | 0.000% [0.000%, 0.000%] | +1.108% |
| 0-10% | 682,860 | 3,147 | 40 | 12,224 | 1.790% [1.317%, 2.387%] | 0 | 0.000% [0.000%, 0.000%] | +1.790% |
| 10-30% | 72,130 | 319 | 25 | 1,336 | 1.852% [1.026%, 3.082%] | 0 | 0.000% [0.000%, 0.000%] | +1.852% |
| >30% | 26,144 | 117 | 11 | 615 | 2.352% [0.859%, 4.937%] | 0 | 0.000% [0.000%, 0.000%] | +2.352% |

---

## 4. Cycle Mensuel en Clay Loam : Contraste GMIA >30 % vs 0 %

| Mois | SMAP >30% [IC 95%] | SMAP 0% [IC 95%] | Contraste SMAP (>30% − 0%) | OL >30% | OL 0% | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 01 | 7.736% [4.172%, 13.591%] | 5.251% [3.519%, 6.854%] | +2.485% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |
| 02 | 5.003% [2.424%, 9.957%] | 3.260% [2.112%, 4.358%] | +1.743% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |
| 03 | 4.271% [2.443%, 8.607%] | 3.866% [2.662%, 5.261%] | +0.405% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |
| 04 | 3.563% [2.110%, 6.949%] | 3.147% [1.867%, 4.091%] | +0.416% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |
| 05 | 3.880% [0.923%, 10.382%] | 2.220% [1.490%, 3.024%] | +1.660% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |
| 06 | 4.264% [0.669%, 11.705%] | 2.486% [1.323%, 3.683%] | +1.778% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 07 | 5.010% [1.194%, 11.652%] | 1.581% [0.663%, 2.707%] | +3.429% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 08 | 4.348% [1.224%, 10.265%] | 2.230% [1.324%, 3.292%] | +2.118% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 09 | 3.884% [0.551%, 9.587%] | 1.847% [0.635%, 3.427%] | +2.038% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 10 | 4.133% [1.459%, 9.200%] | 2.118% [0.896%, 3.307%] | +2.016% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 11 | 3.934% [1.909%, 8.219%] | 5.147% [2.689%, 7.278%] | -1.214% | 0.000% | 0.000% | 252 / 329 | 9 / 15 |
| 12 | 5.298% [2.177%, 10.023%] | 5.404% [2.757%, 7.424%] | -0.106% | 0.000% | 0.000% | 252 / 331 | 9 / 15 |

---

## 5. Analyse de Sensibilité aux Seuils

| Configuration | Seuils appliqués | SMAP JJA >30% | SMAP JJA 0% | Contraste SMAP JJA | OL JJA >30% |
|:---|:---|:---:|:---:|:---:|:---:|
| default | ΔSM > 0.03 m³/m³, Rain < 1.0 mm | 4.531% | 2.103% | +2.428% | 0.000% |
| sensitivity_lax | ΔSM > 0.02 m³/m³, Rain < 1.0 mm | 8.867% | 5.935% | +2.932% | 0.000% |
| sensitivity_strict | ΔSM > 0.04 m³/m³, Rain < 0.5 mm | 2.379% | 0.646% | +1.733% | 0.000% |

---

## 6. Référence aux Figures et Données Produites

- **Figure F6** : `figures/descriptive/F6_dynamic_irrigation_wetting_events.png` (version EN : `figures/descriptive/en/F6_dynamic_irrigation_wetting_events_en.png`)
- **Table complète par strate** : `reports/phase2a_irrigation_diagnostic.csv`
- **Table d'analyse de sensibilité** : `reports/phase2a_irrigation_diagnostic_sensitivity.csv`
