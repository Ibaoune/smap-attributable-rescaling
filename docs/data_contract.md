# Contrat de Données : Table de Modélisation SAR

Ce contrat définit les spécifications structurelles, sémantiques et de types de données
pour les tables de modélisation du projet **SMAP Attributable Rescaling (SAR)**.

## 1. Format et Organisation
- **Format de stockage** : Apache Parquet avec compression Snappy / ZSTD.
- **Partitionnement physique** : Un fichier par année civile :
  `data/processed/modeling_table/sar_modeling_table_{YYYY}.parquet` pour $YYYY \in \{2016, 2017, 2018, 2019, 2020\}$.
- **Unité d'observation (granularité)** : Une ligne par triplet **(pixel, jour, passage)** où une observation SMAP brute recommandée existe sur le domaine NorthMor.

---

## 2. Dictionnaire des Variables

### A. Clés Primaires et Métadonnées Spatio-Temporelles
| Nom de colonne | Type | Unité | Plage / Valeurs | Description |
| :--- | :--- | :---: | :---: | :--- |
| `date` | `date32` / `string` | — | `YYYY-MM-DD` | Date civile UTC de l'observation |
| `year` | `int16` | — | 2016 à 2020 | Année |
| `month` | `int8` | — | 1 à 12 | Mois |
| `day` | `int8` | — | 1 à 31 | Jour du mois |
| `doy` | `int16` | — | 1 à 366 | Jour de l'année (Day of Year) |
| `pass` | `string` | — | `'AM'`, `'PM'` | Passage orbital SMAP (AM ~06h00 locale, PM ~18h00 locale) [Groupe removable] |
| `lat` | `float32` | °N | [30.00, 36.00] | Latitude du centre de maille NorthMor (vérifiée) |
| `lon` | `float32` | °E | [-10.00, -1.00] | Longitude du centre de maille NorthMor (vérifiée) |
| `y_index` | `int16` | indice | 0 à 120 | Indice de ligne sur la grille NorthMor (sud vers nord) |
| `x_index` | `int16` | indice | 0 à 180 | Indice de colonne sur la grille NorthMor (ouest vers est) |
| `spatial_block_id` | `string` | — | Ex. `'33_-7'` | Identifiant du bloc spatial de 1° × 1° |
| `spatial_fold` | `int8` | indice | 0 à 4 | Numéro du pli de validation croisée spatiale (0 à 4) |
| `temporal_split` | `string` | — | `'train'`, `'val'`, `'test'` | Rôle temporel : train (2016–2018), val (2019), test (2020 scellé) |

### B. Cible d'Apprentissage et Composantes Cible (NON-PRÉDICTEURS)
| Nom de colonne | Type | Unité | Statut | Description |
| :--- | :--- | :---: | :---: | :--- |
| `d` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | **CIBLE** | $d = \text{smap\_sm} - \text{ssm\_ol\_daily}$ |
| `smap_sm` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | **NON-PRÉDICTEUR** | Humidité SMAP brute reprojectée (composante de $d$, strictement exclue des prédicteurs) |
| `ssm_ol_daily` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | **NON-PRÉDICTEUR** | Humidité OL journalière au jour $t$ (composante de $d$, strictement exclue des prédicteurs) |

### C. Prédicteurs : Groupe `removable` (À Retirer de SMAP)
| Nom de colonne | Type | Unité | Description |
| :--- | :--- | :---: | :--- |
| `elevation` | `float32` | $\mathrm{m}$ | Altitude SRTM agrégée LDT |
| `slope_deg` | `float32` | degrés | Pente topographique continue calculée à partir de l'altitude |
| `slopetype` | `int8` | catégorie (1–9) | Indice de type de pente LIS (`lis_input_NorthMor_5km.nc`) |
| `soil_texture_fraction_{1..16}` | `float32` | fraction [0, 1] | Fractions texturales STATSGO des 16 classes de sol |
| `land_cover_fraction_{1..20}` | `float32` | fraction [0, 1] | Fractions d'occupation du sol IGBP-MODIS des 20 classes |
| `smap_vegetation_opacity` | `float32` | sans | Épaisseur optique de la végétation (VOD) à bande L |
| `smap_vegetation_water_content` | `float32` | $\mathrm{kg/m}^2$ | Contenu en eau de la végétation (VWC) estimé par SMAP |
| `smap_surface_temperature` | `float32` | $\mathrm{K}$ | Température effective de surface SMAP |
| `smap_roughness_coefficient` | `float32` | sans | Coefficient de rugosité de surface $h$ de l'algorithme tau-omega |
| `smap_static_water_body_fraction` | `float32` | fraction [0, 1] | Fraction de plans d'eau statiques dans la cellule SMAP |
| `smap_surface_flag` | `uint16` / `int32` | bitmask | Drapeaux d'état de surface SMAP |
| `pass` | `string` | `'AM'`, `'PM'` | Passage orbital satellitaire |

### D. Prédicteurs : Groupe `retained` (À Conserver dans l'Observation)
| Nom de colonne | Type | Unité | Description |
| :--- | :--- | :---: | :--- |
| `irrigfrac_gmia` | `float32` | fraction [0, 1] | Fraction d'irrigation GMIA (conservée pour ne pas effacer le signal d'irrigation) |
| `precip_mswep` | `float32` | $\mathrm{mm/d}$ | Précipitation journalière MSWEP v2.8 (référence alternative) |
| `delta_precip_imerg_mswep` | `float32` | $\mathrm{mm/d}$ | Écart forçage pluie : `precip_imerg_daily` − `precip_mswep` |
| `precip_1d_sum` | `float32` | $\mathrm{mm}$ | Pluie antécédente sur 1 jour (strictly $t-1$) |
| `precip_3d_sum` | `float32` | $\mathrm{mm}$ | Pluie antécédente cumulée sur 3 jours (strictly ending at $t-1$) |
| `precip_7d_sum` | `float32` | $\mathrm{mm}$ | Pluie antécédente cumulée sur 7 jours (strictly ending at $t-1$) |
| `days_since_rain_gt_1mm` | `int16` | jours | Nombre de jours consécutifs écoulés sans pluie > 1 mm (strictly ending at $t-1$) |

### E. Prédicteurs : Groupe `context` (Statut Testé par Sensibilité)
| Nom de colonne | Type | Unité | Description |
| :--- | :--- | :---: | :--- |
| `doy_sin` | `float32` | sans | Composante saisonnière annuelle : $\sin(2\pi \cdot \mathrm{DOY} / 365.25)$ |
| `doy_cos` | `float32` | sans | Composante saisonnière annuelle : $\cos(2\pi \cdot \mathrm{DOY} / 365.25)$ |
| `ssm_lag1` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Humidité OL couche 1 au jour $t-1$ |
| `sm_layer2_lag1` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Humidité OL couche 2 au jour $t-1$ |
| `sm_layer3_lag1` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Humidité OL couche 3 au jour $t-1$ |
| `sm_layer4_lag1` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Humidité OL couche 4 au jour $t-1$ |
| `rzsm_lag1` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Humidité moyenne zone racinaire OL au jour $t-1$ |
| `et_lag1` | `float32` | $\mathrm{kg/m}^2/\mathrm{s}$ ou $\mathrm{mm/d}$ | Évapotranspiration totale OL au jour $t-1$ |
| `surface_runoff_lag1` | `float32` | $\mathrm{kg/m}^2/\mathrm{s}$ ou $\mathrm{mm/d}$ | Ruissellement de surface OL au jour $t-1$ |
| `baseflow_lag1` | `float32` | $\mathrm{kg/m}^2/\mathrm{s}$ ou $\mathrm{mm/d}$ | Écoulement de base OL au jour $t-1$ |
| `lai_lag1` | `float32` | $\mathrm{m}^2/\mathrm{m}^2$ | LAI de l'OL au jour $t-1$ |
| `ssm_rel_wilt_ref_lag1` | `float32` | indice [0, 1.5] | Humidité relative $(\text{ssm\_lag1} - \theta_{\text{wlt}}) / (\theta_{\text{ref}} - \theta_{\text{wlt}})$ via SOILPARM.TBL |

### F. Colonnes de Masques Booléens (Contrôles LIS & SMAP)
| Nom de colonne | Type | Description |
| :--- | :---: | :--- |
| `qc_smap_recommended` | `bool` | `True` si bit 0 de `smap_retrieval_qual_flag` == 0 (qualité recommandée SMAP) |
| `qc_lis_rain` | `bool` | `True` si pluie incidente Noah-MP > $3\cdot 10^{-6}\,\mathrm{kg/m}^2/\mathrm{s}$ |
| `qc_lis_frozen` | `bool` | `True` si sol gelé dans Noah-MP ($T_{soil} \le 273.15\,\mathrm{K}$ ou $|smc - sh2o| > 10^{-4}$) |
| `qc_lis_forest` | `bool` | `True` si classe dominante de végétation est une forêt ($\le 4$) |
| `qc_lis_snow` | `bool` | `True` si neige présente ($SWE > 1\,\mathrm{mm}$ ou $SCA > 10^{-4}$) |
| `qc_lis_high_veg` | `bool` | `True` si fraction de végétation verte $shdfac > 0,9$ |
| `qc_lis_sat_limit` | `bool` | `True` si observation trop proche de la saturation : $\theta_{sat} - \text{smap\_sm} < 0,02$ |
| `qc_lis_wilt_limit` | `bool` | `True` si observation trop proche du flétrissement : $\text{smap\_sm} - \theta_{wlt} < 0,02$ |
| `qc_lis_vwc_gt_5` | `bool` | `True` si contenu en eau de la végétation SMAP $VWC > 5,0\,\mathrm{kg/m}^2$ |

### G. Diagnostics EnKF de Référence (Copie pour Intercomparaison et Figures)
| Nom de colonne | Type | Unité | Description |
| :--- | :--- | :---: | :--- |
| `daily_innovation_nocdf` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Innovation journalière EnKF sans CDF matching |
| `daily_signed_increment_nocdf` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Incrément d'analyse signé EnKF sans CDF matching |
| `daily_innovation_cdf_v2` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Innovation journalière EnKF avec CDF matching v2 |
| `daily_signed_increment_cdf_v2` | `float32` | $\mathrm{m}^3/\mathrm{m}^3$ | Incrément d'analyse signé EnKF avec CDF matching v2 |
