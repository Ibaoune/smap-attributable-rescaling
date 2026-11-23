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

---

## Décision 9 : Diagnostic du Signal d'Irrigation Dynamique (Événements d'Humidification sans Pluie)
- **Contexte** : Dans la Phase 2a, le contraste statique GMIA >30% vs 0% sur l'écart $d = \text{SMAP} - \text{OL}$ présente un intervalle de confiance recouvrant zéro en Clay Loam. Pour quantifier la dynamique événementielle des hausses d'humidité indépendamment du biais statique moyen, un diagnostic temporel d'« humidification sans pluie » est appliqué.
- **Décision** :
  1. **Définition de l'événement** : Paires d'observations satellitaires SMAP consécutives du même passage orbital (AM–AM ou PM–PM), séparées d'au plus 3 jours ($1 \le \Delta t \le 3$), satisfaisant simultanément :
     - Hausse de l'humidité du sol : $\Delta\text{SM} = \text{SM}(t) - \text{SM}(t - \Delta t) > \text{seuil}_{\text{SM}}$ (par défaut 0,03 m³/m³).
     - Précipitation IMERG cumulée sur l'intervalle $[t - \Delta t, t]$ : $P_{\text{cum}} < \text{seuil}_{P}$ (par défaut 1,0 mm).
  2. **Témoin modèle (Open Loop Noah-MP)** : Le même diagnostic est appliqué strictement sur la série temporelle de l'OL (`ssm_ol_daily`). Noah-MP étant simulé sans module d'irrigation (`noirr`), ce témoin permet de valider le comportement hydrologique attendu (zéro réhumidification sans apport pluviométrique).
  3. **Stratification et incertitudes** : Calcul de la fréquence d'occurrence par classe GMIA (0%, 0-10%, 10-30%, >30%) $\times$ texture $\times$ mois. Intervalles de confiance à 95% calculés par bootstrap par blocs spatiaux de 1° (500 itérations), avec décompte systématique des mailles et des blocs distincts.
  4. **Sensibilité** : Test de deux autres jeux de seuils (laxe : 0,02 m³/m³ et 1,0 mm ; strict : 0,04 m³/m³ et 0,5 mm) configurés dans `configs/irrigation_diagnostic.yaml`.
- **Statut** : Acté.

---

## Décision 10 : Protocole d'Entraînement et de Rescaling Attribuable Phase 2b (Random Forest)
- **Contexte** : Apprentissage du premier modèle de correction supervisée explicable alternative au CDF matching (Jalon H4).
- **Décision** :
  1. **Population d'entraînement (2016–2018)** : Observations passant l'ensemble des contrôles de qualité (QC) de LIS DA à l'exclusion de la règle de flétrissement.
  2. **Sous-échantillonnage stratifié spatialement** : Pour pallier l'hétérogénéité orbitale et la disparité des masquages QC entre cellules, un échantillonnage aléatoire stratifié uniforme de 20 observations par pixel est appliqué sur 2016–2018 (10 par an maximum).
  3. **Validation croisée spatiale équilibrée** : 5 plis spatiaux constitués de blocs géographiques de 1° regroupés par compacité spatiale et équilibrés en effectif de mailles.
  4. **Formulation du Rescaling Attribuable (SAR)** :
     - Modèle principal (SAR main) : $\text{SMAP}^* = \text{SMAP} - (\Phi_{\text{removable}} + \phi_0)$. Le groupe `retained` (dont `irrigfrac_gmia` et les écarts pluviométriques) est conservé intact dans l'observation corrigée.
     - Variante de sensibilité : $\text{SMAP}^*_{\text{ctx}} = \text{SMAP} - (\Phi_{\text{removable}} + \Phi_{\text{context}} + \phi_0)$, permettant de quantifier l'impact du retrait du contexte temporel/modèle.
  5. **Recalcul du flétrissement et protocole d'évaluation** : La règle $\text{SMAP}^* < \theta_{\text{wlt}} + 0.02$ est recalculée pour chaque observation corrigée sur l'année de validation 2019 (`split == 'val'`). Évaluation strictement identique au protocole Phase 2a (échantillons spécifique et commun à B0/B1/B2/RF/RF-ctx). Année 2020 scellée.
- **Statut** : Acté.

---

## Décision 11 : Remplacement du Témoin OL par Différence de Différences (DiD) et Tendance Continue
- **Contexte** : Dans la simulation Noah-MP Open Loop, l'irrigation n'est pas modélisée (`noirr`). Le témoin OL d'humidification sans pluie est donc nul par construction (0,000 % d'événements) et n'apporte aucune information empirique supplémentaire.
- **Décision** :
  1. **Différence de différences (DiD)** : Pour contrôler les effets fixes spatiaux de texture et les variations saisonnières climatiques, le signal dynamique est quantifié par :
     $$\text{DiD} = [f(>30\,\%) - f(0\,\%)]_{\text{JJA}} - [f(>30\,\%) - f(0\,\%)]_{\text{DJF}}$$
     à texture fixée, avec intervalles de confiance à 95 % par bootstrap sur blocs de 1°.
  2. **Test de tendance continue** : Régression linéaire de la fréquence d'humidification sans pluie en fonction de la fraction GMIA continue ($[0, 1]$), par texture et par saison, rapportant la pente $\beta$ (%/fraction GMIA) et la corrélation de Pearson $r$, avec intervalles de confiance à 95 % par bootstrap sur blocs de 1°.
- **Statut** : Acté et intégré dans `scripts/08_diagnose_dynamic_irrigation.py`.

---

## Décision 12 : Passage à LightGBM et Attribution TreeSHAP Exacte (Phase 2c)
- **Contexte** : Le Random Forest avec TreeExplainer Python et approximations de Saabas limitait la taille d'entraînement (20 obs/pixel) et présentait des déviations d'additivité.
- **Décision** :
  1. **Algorithme principal** : LightGBM (`LGBMRegressor`, 300 estimateurs, `learning_rate=0.05`, `num_leaves=63`). Le Random Forest est conservé uniquement à titre de comparaison historique.
  2. **Taille d'échantillon d'entraînement** : Échantillon stratifié de 100 observations par pixel tirées uniformément parmi les observations non-wilt QC-valides sur 2016–2018 (soit 1 425 000 observations sur 14 250 pixels).
  3. **Attribution TreeSHAP exacte** : Calcul natif C++ par `model.booster_.predict(X, pred_contrib=True)` sur l'intégralité des 3 660 799 observations de l'année de validation 2019. Aucune approximation de Saabas n'est autorisée. L'additivité exacte de Lundberg ($\sum \phi_j + \phi_0 = \hat{y}$) est vérifiée numériquement (déviation maximale $\le 10^{-15}$).
- **Statut** : Acté et intégré dans `scripts/10_train_phase2c_lightgbm.py`.

---

## Décision 13 : Partitionnement du Contexte et Variantes de Rescaling (V1, V2, V3)
- **Contexte** : Dans `configs/predictors.yaml`, le groupe contexte regroupe à la fois des représentations temporelles/statiques d'état et des flux dynamiques de processus.
- **Décision** :
  1. **Scission du groupe contexte** :
     - `context_representation` (4 variables) : `doy_sin`, `doy_cos`, `ssm_lag1`, `ssm_rel_wilt_ref_lag1`.
     - `context_process` (8 variables) : `sm_layer2_lag1`, `sm_layer3_lag1`, `sm_layer4_lag1`, `rzsm_lag1`, `et_lag1`, `surface_runoff_lag1`, `baseflow_lag1`, `lai_lag1`.
  2. **Trois variantes de correction additive** (toutes avec $\phi_0$) :
     - **V1 (removable)** : $\text{SMAP}^*_{\text{V1}} = \text{SMAP} - (\Phi_{\text{removable}} + \phi_0)$.
     - **V2 (removable + context_representation)** : $\text{SMAP}^*_{\text{V2}} = \text{SMAP} - (\Phi_{\text{removable}} + \Phi_{\text{context\_rep}} + \phi_0)$.
     - **V3 (removable + all context)** : $\text{SMAP}^*_{\text{V3}} = \text{SMAP} - (\Phi_{\text{removable}} + \Phi_{\text{context\_rep}} + \Phi_{\text{context\_proc}} + \phi_0)$ (analogue de RF-ctx).
  3. **Métriques requises** : Pour chaque variante, rapport systématique de la correction moyenne appliquée $\overline{\Delta\text{SMAP}}$ et du biais résiduel $\overline{d_{\text{res}}}$ par saison (DJF, MAM, JJA, SON) et en annuel, en plus du protocole complet 2a.
- **Statut** : Acté et intégré dans `configs/predictors.yaml` et `scripts/10_train_phase2c_lightgbm.py`.

---

## Décision 14 : Reconstruction LDT de la Référence CDF B2-train (2016-01-01 à 2019-01-01)
- **Contexte** : La référence historique B2 (2016–2021) englobait la période d'évaluation 2019 et la période scellée 2020.
- **Décision** : Reconstruire avec LDT les fichiers CDF observation et modèle strictement sur la période d'entraînement 2016-01-01 00:00:00 à 2019-01-01 00:00:00 dans `experiments/ldt_b2_train/` (accord utilisateur explicite reçu). Cette nouvelle référence, notée **B2-train**, élimine toute fuite temporelle et est intégrée à l'ensemble des comparaisons du protocole 2a.
- **Statut** : Exécuté et validé.
