# Scope du Projet : SMAP Attributable Rescaling (SAR)

## 1. Contexte Scientifique
L'assimilation de données satellitaires d'humidité du sol (SMAP) dans les modèles de surface (LSM, ici Noah-MP sous LIS) requiert classiquement une étape de mise à l'échelle (rescaling) pour éliminer les biais systématiques entre l'observation et la climatologie du modèle.
L'approche standard par **CDF matching** (appariement de fonctions de distribution cumulative) présente plusieurs limites documentées :
1. Elle est purement empirique et statique, traitant le biais comme une distorsion stationnaire quantile par quantile.
2. Elle ignore les mécanismes physiques à l'origine du biais (texture du sol, type de végétation, régimes de précipitations, dynamique d'assèchement, irrigation non simulée).
3. Elle efface aveuglément les signaux physiques non simulés par le modèle, comme l'irrigation agricole.
4. Elle introduit des distorsions numériques aux queues de distribution (flétrissement ou saturation).

Le projet **SMAP Attributable Rescaling (SAR)** a pour ambition de concevoir, entraîner et évaluer une correction de biais **explicable, conditionnelle et attribuable**, capable de prédire l'écart SMAP–modèle à partir de covariables physiques réparties en composantes à retirer, à conserver et contextuelles.

---

## 2. Domaine Géographique et Grille
- **Nom du domaine** : NorthMor (Nord du Maroc).
- **Dimensions de grille** : 121 rangées de latitude × 181 colonnes de longitude.
- **Résolution spatiale** : 0,05° (environ 5 km).
- **Limites géographiques vérifiées** (centres de mailles issus de `lis_input_NorthMor_5km.nc`) :
  - **Latitude** : 30,00°N à 36,00°N (coin SW : 30,0°N ; coin NE : 36,0°N)
  - **Longitude** : -10,00°E à -1,00°E (coin SW : -10,0°E ; coin NE : -1,0°E)
- **Support terrestre actif** : 16 030 mailles continentales valides (masque terre strict de NorthMor).

---

## 3. Période Temporelle
- **Horizon complet** : 2016-01-01 au 2020-12-31 (5 années complètes, soit 1 827 jours civils).
- **Disponibilité SMAP** : 1 792 jours disponibles sur le HPC (35 jours manquants documentés, dont 33 jours lors de l'anomalie instrumentale de juin-juillet 2019).

---

## 4. Cible d'Apprentissage
- **Variable cible principale** : Écart signé $d = \text{SMAP}_{\text{brut}} - \text{SSM}_{\text{OL daily}}$ ($\mathrm{m}^3/\mathrm{m}^3$).
- **Couche de sol** : Couche superficielle 1 (0–10 cm).
- **Passages satellitaires** : Séparation stricte des orbites AM (descendante, ~06:00 heure locale) et PM (ascendante, ~18:00 heure locale).
- **Agrégation temporelle de l'OL** : Moyenne journalière 24h (`SoilMoist_tavg`), l'OL horaire étant reporté (voir `docs/decisions.md`).

---

## 5. Groupes de Prédicteurs (Attribution Fonctionnelle)
1. **Groupe `removable` (retiré de SMAP)** :
   Prédicteurs associés aux artefacts de mesure satellitaire ou aux biais physiographiques systématiques :
   - Altitude SRTM, pente topographique continue, type de pente LIS (`SLOPETYPE`).
   - Fractions texturales du sol STATSGO (16 classes continues).
   - Fractions d'occupation du sol MODIS-IGBP (20 classes continues).
   - Variables de restitution SMAP : VOD (`vegetation_opacity`), VWC (`vegetation_water_content`), température de surface (`surface_temperature`), rugosité (`roughness_coefficient`), fraction d'eau libre (`static_water_body_fraction`), flags de surface (`surface_flag`).
   - Orbite du passage : `pass` (AM / PM).
2. **Groupe `retained` (conservé dans l'observation)** :
   Prédicteurs porteurs de signaux physiques réels non modélisés par l'OL ou d'erreurs de forçage :
   - Fraction d'irrigation GMIA (`irrigfrac_gmia`) : conservée pour ne pas effacer le signal d'irrigation comme le fait le CDF matching.
   - Précipitation MSWEP v2.8 (`precip_mswep`) et écart de forçage `delta_precip_imerg_mswep`.
   - Cumuls de pluie antécédents à 1, 3 et 7 jours (`precip_1d_sum`, `precip_3d_sum`, `precip_7d_sum`).
   - Nombre de jours consécutifs sans pluie (`days_since_rain_gt_1mm`).
3. **Groupe `context` (statut testé par sensibilité)** :
   Variables d'état et de contexte temporel dont l'impact sera quantifié :
   - Harmoniques saisonnières : `doy_sin` et `doy_cos`.
   - États et flux antécédents de l'Open Loop à $t-1$ : `ssm_lag1`, `sm_layer2_lag1` à `sm_layer4_lag1`, `rzsm_lag1`, `et_lag1`, `surface_runoff_lag1`, `baseflow_lag1`, `lai_lag1`.
   - Indice d'humidité relative par rapport au flétrissement et à la capacité au champ à $t-1$ (`ssm_rel_wilt_ref_lag1`).

---

## 6. Protocole Anti-Fuite et Découpage
- **Principe d'étanchéité temporelle** : Aucun prédicteur ne contient l'état de l'OL du jour $t$ ni l'observation `smap_sm`. Tous les états OL et indices dérivés sont strictement indexés à $t-1$.
- **Partition temporelle** : Entraînement et validation croisée sur 2016–2019 ; scellement strict de l'année 2020 comme jeu de test indépendant final.
- **Partition spatiale** : 5 plis spatiaux réguliers équilibrés par blocs de 1° × 1°, interdisant tout découpage aléatoire pixel-jour.
