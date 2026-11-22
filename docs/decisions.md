# Registre des Décisions d'Architecture et Méthodologiques (ADR)

Ce document consigne les décisions méthodologiques formelles actées pour le projet
**SMAP Attributable Rescaling (SAR)**, ainsi que les points ouverts reportés.

---

## Décision 1 : Cible d'apprentissage signée ($d = \text{SMAP}_{\text{brut}} - \text{SSM}_{\text{OL}}$)
- **Contexte** : Le CDF matching produit une variable transformée sans référence directe à la physique du modèle. Pour rendre la correction attribuable et interprétable, nous prédisons directement l'écart de biais.
- **Décision** : La cible est la différence algébrique signée :
  $$d = \text{SMAP}_{\text{brut}} - \text{SSM}_{\text{OL}}$$
  exprimée en $\mathrm{m}^3/\mathrm{m}^3$.
- **Rationnel** : Un signe positif ($d > 0$) indique que le radiomètre SMAP voit un sol plus humide que l'Open Loop ; un signe négatif ($d < 0$) indique un sol vu plus sec. Ce formalisme permet de décomposer l'erreur en composantes attribuables (ex. forçage de pluie sous-estimé, végétation dense, irrigation omise).
- **Statut** : Acté.

---

## Décision 2 : Utilisation de la moyenne journalière de l'Open Loop (`SoilMoist_tavg`)
- **Contexte** : Dans la configuration `OPL_noirr_2016_2020` de Noah-MP sous LIS, les sorties n'ont été écrites qu'au pas journalier (`1da`, moyenne sur 24 heures). Les sorties horaires permettant d'isoler l'état du modèle à 06h00 (passage AM) et 18h00 (passage PM) ne sont disponibles sur le HPC que pour l'année 2018 (`sensitivity_2018/OL_ctrl_hourly_2018`).
- **Décision** : Pour couvrir l'ensemble des 5 années (2016–2020), nous utilisons la moyenne journalière 24h `SoilMoist_tavg` (couche 1, 0–10 cm) comme référence de l'OL.
- **POINT OUVERT (Reporté)** : L'extraction de l'OL à l'heure exacte du passage satellitaire (06h / 18h) est reportée. L'audit réalisé sur l'année 2018 (`smap-xai-audit/reports/audits/diag2018_e6_clean/ol_smap_overpass_sampling_bias.csv`) démontre que le biais d'échantillonnage diurne moyen de l'OL est faible :
  - 06h00 vs moyenne 24h : $+0,0001$ à $+0,0034\,\mathrm{m}^3/\mathrm{m}^3$ selon les saisons et textures.
  - 18h00 vs moyenne 24h : $-0,0002$ à $-0,0036\,\mathrm{m}^3/\mathrm{m}^3$.
  - Amplitude diurne globale : $0,0005$ à $0,0070\,\mathrm{m}^3/\mathrm{m}^3$.
  Ce point pourra être réexaminé si une relance complète de l'OL horaire sur 5 ans est effectuée ultérieurement.
- **Statut** : Acté (moyenne journalière retenue) avec point ouvert documenté.

---

## Décision 3 : Contrôle Qualité (QC) et Filtrage des Observations
- **Contexte** : LIS applique deux niveaux de filtrage :
  1. Le plugin d'observation SMAP (`read_NASASMAPsm.F90`) filtre sur le bit 0 du `retrieval_qual_flag == 0` et élimine les pixels où VWC > 5 kg/m².
  2. Le module de pré-assimilation Noah-MP (`noahmp401_qc_soilmobs.F90`) rejette les observations en cas de pluie incidente ($> 3\cdot 10^{-6}\,\mathrm{kg/m}^2/\mathrm{s}$), sol gelé, neige, forêt (type IGBP $\le 4$), végétation verte dense (`shdfac` $> 0,9$), ou si l'observation est trop proche des bornes physiques du sol ($\text{smobs} - \text{smcwlt} < 0,02$ ou $\text{smcmax} - \text{smobs} < 0,02$).
- **Décision** :
  - Pour définir la validité de la cible SMAP, **seul le flag de qualité SMAP** (bit 0 de `retrieval_qual_flag == 0`, qualité recommandée) est utilisé comme critère d'inclusion.
  - Les contrôles environnementaux de LIS (`noahmp401_qc_soilmobs.F90`) **ne censurent pas** l'échantillon d'apprentissage de SAR.
  - Au lieu de cela, chaque contrôle LIS est reproduit rigoureusement et conservé sous forme de **colonne de masque booléen** (`qc_lis_rain`, `qc_lis_frozen`, `qc_lis_forest`, `qc_lis_snow`, `qc_lis_high_veg`, `qc_lis_sat_limit`, `qc_lis_wilt_limit`, `qc_lis_vwc_gt_5`).
- **Rationnel** : Cela permet aux modèles SAR d'apprendre sur la distribution complète des observations satellitaires recommandées, tout en rendant possible l'audit précis de la part des données rejetées par LIS et l'analyse de sensibilité.
- **Statut** : Acté.

---

## Décision 4 : Reprojection spatiale par plus proche voisin identique à LIS
- **Contexte** : Le système LIS assimile les données SMAP Enhanced (EASE-Grid 2.0 9 km) en projetant chaque maille du domaine cible (NorthMor 121 × 181, 0,05°) vers la cellule SMAP la plus proche via la sous-routine Fortran `neighbor_interp` (`neighbor_interp_input_withgrid.F90`).
- **Décision** : Le module `src/sar/regrid.py` reproduit de manière déterministe et exacte l'arrondi au plus proche voisin de LIS (`nint(xi)`, `nint(yi)` dans le repère de grille).
- **Validation** : Cette implémentation doit être validée empiriquement par comparaison directe avec les fichiers `DAOBS` issus de la simulation `DA_nocdf_noirr_2016_2020` (écart attendu strictement nul, aux arrondis de troncature près).
- **Statut** : Acté.

---

## Décision 5 : Séparation stricte des passages AM et PM
- **Contexte** : Dans `read_NASASMAPsm.F90` (lignes 280–298), LIS DA applique une règle asymétrique : le passage AM (06h) est prioritaire, et le passage PM (18h) n'est injecté que là où l'observation AM est absente.
- **Décision** : Pour le projet SAR, les deux passages orbitaux sont traités de manière autonome et symétrique. Chaque pixel-jour-passage constitue une ligne distincte dans la table de modélisation avec un attribut `pass` $\in \{\text{'AM'}, \text{'PM'}\}$.
- **Rationnel** : La physique de surface (température, gradients d'humidité, rosée matinale, stratification de la couche limite) diffère profondément entre 06h00 et 18h00. Écraser un passage au profit de l'autre biaiserait l'apprentissage de la correction.
- **Statut** : Acté.

---

## Décision 6 : Structuration des Prédicteurs en Trois Groupes Fonctionnels
- **Contexte** : Le rescaling par CDF matching classique projette indistinctement toutes les variations de l'observation sur la climatologie du modèle, effaçant ainsi des signaux réels comme l'irrigation agricole non simulée.
- **Décision** : Les prédicteurs de SAR sont explicitement structurés en trois groupes :
  1. `removable` (à retirer de SMAP) : Variables liées aux artefacts instrumentaux et biais systématiques (altitude, pente continue, slopetype, fractions de texture, fractions d'occupation du sol, VOD, VWC, rugosité, Tsfc, fraction d'eau libre, surface flag, passage AM/PM).
  2. `retained` (à conserver dans l'observation) : Signaux physiques réels absents du modèle et incertitudes de forçage :
     - `irrigfrac_gmia` : **conservé** impérativement pour préserver l'incrément d'irrigation détecté par SMAP et empêcher son écrasement.
     - `precip_mswep`, `delta_precip_imerg_mswep` : écarts de forçage pluviométrique.
     - `precip_1d_sum`, `precip_3d_sum`, `precip_7d_sum`, `days_since_rain_gt_1mm` : mémoire hydrique antécédente.
  3. `context` (statut testé par analyse de sensibilité) : Variables temporelles et d'état antécédent du modèle (`doy_sin`, `doy_cos`, états OL à $t-1$, `ssm_rel_wilt_ref` à $t-1$). Leur statut définitif (conservé ou retiré) sera tranché expérimentalement.
- **Statut** : Acté.

---

## Décision 7 : Règle Stricte Anti-Fuite (Étanchéité Temporelle et Non-Prédicteurs)
- **Contexte** : Dans un cadre supervisé visant à prédire $d = \text{SMAP}_{\text{brut}} - \text{SSM}_{\text{OL daily}}$, l'utilisation directe de grandeurs contemporaines au jour $t$ issues du modèle ou de la cible crée une fuite d'information directe ou un biais d'anticipation.
- **Décision** :
  1. **États du modèle à $t-1$** : Aucun prédicteur ne contient l'état de l'OL du jour $t$. Toutes les variables d'état et flux de l'OL (`ssm_lag1`, `sm_layer2_lag1` à `sm_layer4_lag1`, `rzsm_lag1`, `et_lag1`, `surface_runoff_lag1`, `baseflow_lag1`, `lai_lag1`, `ssm_rel_wilt_ref_lag1`) sont décalées d'au moins 1 jour ($t-1$).
  2. **Interdiction de `smap_sm` comme prédicteur** : L'observation SMAP brute du jour $t$ intervient exclusivement comme composante de la cible $d$. Elle n'est en aucun cas un prédicteur explicatif.
  3. **Précipitations du jour $t$** : La précipitation totale journalière IMERG du jour $t$ n'est pas connue à 06h00 (passage AM). Par conséquent, les cumuls de pluie utilisés comme prédicteurs sont strictement antécédents (`precip_1d_sum`, `precip_3d_sum`, `precip_7d_sum` arrêtés à $t-1$). Si la pluie du jour $t$ devait être utilisée, elle exigerait d'être tronquée à l'heure du passage satellitaire via les forçages demi-horaires IMERG (05h00 UTC pour AM, 17h00 UTC pour PM) ; à défaut, elle est strictement exclue des prédicteurs pour éliminer toute fuite.
  4. **Composantes de la cible dans la table** : `ssm_ol_daily` et `smap_sm` sont stockées dans la table Parquet à des fins de calcul de diagnostics et d'audit, mais sont formellement marquées « non-prédicteur » dans `docs/data_contract.md`.
- **Statut** : Acté.

---

## Décision 8 : Source autoritaire pour la fraction irriguée (FAO GMIA v5)
- **Contexte** : Dans le fichier maître de paramètres d'entrée de LIS `NoahMP_Morocco/data/lis_input/lis_input_NorthMor_5km.nc`, la variable `IRRIGFRAC` est strictement nulle sur l'ensemble des 21 901 mailles du domaine NorthMor (`min=0.0`, `max=0.0`, `mean=0.0`), de même que dans `lis_input_NorthMor_5km_irr.nc`. Cette vacuité s'explique par le fait que les simulations historiques Noah-MP ont été exécutées en configuration non irriguée (`noirr`), sans initialisation de la fraction d'irrigation dans `lis_input`.
- **Décision** : La source faisant autorité pour quantifier la fraction équipée pour l'irrigation est la base de données globale de la FAO : **GMIA v5** (*Global Map of Irrigation Areas*, version 5.0), située dans `NoahMP_Morocco/data/land_params/GMIA/gmia_v5_aei_pct.asc` (matrice binaire optimisée `gmia_v5_aei_pct.bin`).
- **Méthode** :
  - Sous-ensemble extrait sur la fenêtre géographique du domaine NorthMor ($[29.8^\circ\text{N}, 36.2^\circ\text{N}] \times [-10.2^\circ\text{E}, -0.8^\circ\text{E}]$).
  - Interpolation spatiale sur la grille cible 0,05° (`lat`, `lon`) par interpolation bi-linéaire sur grille régulière (`scipy.interpolate.RegularGridInterpolator`, `bounds_error=False`, `fill_value=0.0`).
  - Conversion d'unité : passage du pourcentage $[0, 100\,\%]$ à une fraction adimensionnelle $[0, 1]$ (`irrigfrac_gmia = gmia_5km / 100.0`).
- **Validation** : La grille interpolée présente une moyenne de 0,0255 (2,55 %) sur NorthMor, avec un maximum de 0,9550 (95,50 % dans le bassin du Gharb/Sebou).
- **Statut** : Acté.
