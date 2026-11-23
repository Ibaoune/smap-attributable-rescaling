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

---

## 2. Statistiques Globales (Configuration par Défaut)

- **Nombre total de paires de passages successifs** : 13,878,569
- **Événements détectés sur SMAP (ensemble du domaine et de l'année)** : 261,486 (1.884%, IC 95% : [1.465%, 2.288%])

---

## 3. Différence de Différences (DiD) : [Fréquence >30 % − 0 %] en Été (JJA) vs Hiver (DJF)

Contraste $\Delta_{	ext{DiD}} = [f(>30\%) - f(0\%)]_{	ext{JJA}} - [f(>30\%) - f(0\%)]_{	ext{DJF}}$ avec intervalles de confiance à 95% par bootstrap par blocs de 1° (500 itérations) :

| Texture | Contraste JJA (>30% − 0%) | Contraste DJF (>30% − 0%) | Différence de Différences (DiD) | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---|:---:|:---:|:---:|:---:|:---:|
| Clay Loam | +2.428% [-1.495%, +10.731%] | +1.359% [-2.135%, +6.757%] | **+1.069%** [-2.347%, +5.379%] | 252/331 px | 9/15 blk |
| Loam | +1.244% [-0.357%, +3.905%] | +1.514% [-0.158%, +3.172%] | **-0.270%** [-1.254%, +1.952%] | 117/7467 px | 11/46 blk |
| Sandy Clay Loam | +nan% [+nan%, +nan%] | +nan% [+nan%, +nan%] | **+nan%** [+nan%, +nan%] | 0/613 px | 0/8 blk |
| Sandy Loam | +9.138% [-2.241%, +21.604%] | +5.874% [-1.452%, +14.525%] | **+3.264%** [-0.981%, +7.183%] | 14/312 px | 3/12 blk |
| Clay | +5.127% [-1.858%, +12.717%] | -0.288% [-4.691%, +4.897%] | **+5.414%** [+2.114%, +8.523%] | 26/136 px | 3/7 blk |
| All | +3.095% [+0.360%, +7.595%] | +3.355% [+0.917%, +6.266%] | **-0.260%** [-2.411%, +2.328%] | 433/9412 px | 14/48 blk |

---

## 4. Test de Tendance Continue : Fréquence d'Humidification vs Fraction GMIA Continue

Pente de régression linéaire $eta$ (% de fréquence par 100% de GMIA) et corrélation de Pearson $r$ avec intervalles de confiance à 95% par bootstrap par blocs de 1° :

| Texture | Saison | Paires analysées | Pente linéaire $eta$ [IC 95%] | Corrélation $r$ [IC 95%] | Échantillon spatial |
|:---|:---:|:---:|:---:|:---:|:---:|
| Clay Loam | JJA | 267,286 | +3.277% [-1.699%, +16.027%] | +0.0527 [-0.0306, +0.1767] | 1,056 px, 19 blk |
| Clay Loam | DJF | 273,734 | +3.057% [-2.564%, +11.734%] | +0.0395 [-0.0332, +0.1291] | 1,058 px, 19 blk |
| Clay Loam | All | 1,100,145 | +2.560% [-1.977%, +12.270%] | +0.0384 [-0.0298, +0.1351] | 1,058 px, 19 blk |
| Loam | JJA | 2,577,955 | +3.852% [+0.397%, +10.288%] | +0.0193 [+0.0018, +0.0446] | 11,050 px, 48 blk |
| Loam | DJF | 2,585,372 | +4.244% [+0.452%, +7.789%] | +0.0177 [+0.0018, +0.0350] | 11,050 px, 48 blk |
| Loam | All | 10,510,564 | +4.135% [+0.835%, +8.694%] | +0.0193 [+0.0034, +0.0391] | 11,051 px, 48 blk |
| All | JJA | 3,399,963 | +5.318% [+1.089%, +12.810%] | +0.0489 [+0.0105, +0.1002] | 14,222 px, 50 blk |
| All | DJF | 3,409,951 | +6.133% [+1.794%, +12.093%] | +0.0435 [+0.0105, +0.0803] | 14,249 px, 50 blk |
| All | All | 13,878,569 | +5.176% [+1.055%, +11.692%] | +0.0428 [+0.0088, +0.0819] | 14,250 px, 50 blk |

---

## 5. Cycle Mensuel en Clay Loam : Contraste GMIA >30 % vs 0 %

| Mois | SMAP >30% [IC 95%] | SMAP 0% [IC 95%] | Contraste SMAP (>30% − 0%) | Pixels (>30% / 0%) | Blocs 1° (>30% / 0%) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 01 | 7.736% [3.885%, 13.367%] | 5.251% [3.501%, 7.021%] | +2.485% | 252/331 | 9/15 |
| 02 | 5.003% [2.293%, 9.321%] | 3.260% [2.031%, 4.245%] | +1.743% | 252/331 | 9/15 |
| 03 | 4.271% [2.445%, 8.728%] | 3.866% [2.512%, 5.093%] | +0.405% | 252/331 | 9/15 |
| 04 | 3.563% [2.148%, 6.875%] | 3.147% [1.878%, 4.115%] | +0.416% | 252/331 | 9/15 |
| 05 | 3.880% [1.122%, 10.728%] | 2.220% [1.446%, 2.993%] | +1.660% | 252/331 | 9/15 |
| 06 | 4.264% [0.671%, 13.490%] | 2.486% [1.307%, 3.715%] | +1.778% | 252/329 | 9/15 |
| 07 | 5.010% [0.881%, 14.189%] | 1.581% [0.627%, 2.627%] | +3.429% | 252/329 | 9/15 |
| 08 | 4.348% [1.297%, 9.888%] | 2.230% [1.279%, 3.390%] | +2.118% | 252/329 | 9/15 |
| 09 | 3.884% [0.591%, 10.772%] | 1.847% [0.540%, 3.382%] | +2.038% | 252/329 | 9/15 |
| 10 | 4.133% [1.204%, 10.102%] | 2.118% [0.908%, 3.435%] | +2.016% | 252/329 | 9/15 |
| 11 | 3.934% [1.999%, 8.274%] | 5.147% [2.896%, 7.415%] | -1.214% | 252/329 | 9/15 |
| 12 | 5.298% [2.036%, 10.832%] | 5.404% [2.941%, 7.552%] | -0.106% | 252/331 | 9/15 |

---

## 6. Analyse de Sensibilité aux Seuils

| Configuration | Seuils appliqués | Clay Loam >30% | Clay Loam 0% | Contraste Clay Loam | Loam >30% | Loam 0% | Contraste Loam |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| default | ΔSM > 0.03 m³/m³, Rain < 1.0 mm | 4.531% | 2.103% | +2.428% | 2.352% | 1.108% | +1.244% |
| sensitivity_lax | ΔSM > 0.02 m³/m³, Rain < 1.0 mm | 8.867% | 5.935% | +2.932% | 5.653% | 2.862% | +2.791% |
| sensitivity_strict | ΔSM > 0.04 m³/m³, Rain < 0.5 mm | 2.379% | 0.646% | +1.733% | 1.021% | 0.338% | +0.683% |

---

## 7. Fichiers et Données Produites

- **Figure F6** : `figures/descriptive/F6_dynamic_irrigation_wetting_events.png` (version EN : `figures/descriptive/en/F6_dynamic_irrigation_wetting_events_en.png`)
- **Table complète par strate** : `reports/phase2a_irrigation_diagnostic.csv`
- **Table Différence de Différences (DiD)** : `reports/phase2a_irrigation_did.csv`
- **Table test de tendance continue GMIA** : `reports/phase2a_irrigation_trend.csv`
- **Table d'analyse de sensibilité** : `reports/phase2a_irrigation_diagnostic_sensitivity.csv`
