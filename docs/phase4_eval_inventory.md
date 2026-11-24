# Inventaire des Ressources d'Évaluation (Phase 1 et Pilote DA-V2)

Date de référence : 24 novembre 2026  
Auteur : M. El Aabaribaoune (@um6p)  
Statut : Document de référence en lecture seule pour la Phase 4 (Pilote DA 2019 et extension 2016–2020)

---

## 1. Périmètre et Objectif de l'Inventaire

Ce document dresse l'inventaire exhaustif et vérifié des données d'observation, configurations, modèles hydrodynamiques et jeux de données auxiliaires existant sur le cluster Toubkal pour évaluer :
1. Les simulations historiques de la Phase 1 (OPL, DA-NoCDF, DA-CDF sur 2016–2020) ;
2. Les trois simulations pilotes de l'année 2019 de la Phase 4 (`DA_nocdf_2019`, `DA_cdftrain_2019`, `DA_v2_2019`).

L'inventaire couvre quatre volets spécifiques :
- Les stations hydrométriques de débit (Sebou et Bouregreg) ;
- La configuration et les scripts du routage hydrodynamique HyMAP2 de P1 ;
- Les produits d'évapotranspiration (ET) indépendants disponibles ;
- Les mesures in situ d'humidité du sol et leurs produits de substitution.

---

## 2. Stations de Débit In Situ (Sebou et Bouregreg)

### 2.1 Emplacement et Référentiel sur le Cluster
- **Répertoire principal nettoyé :**  
  `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/streamflow_clean/`
- **Fichiers clés :**
  - `stations.csv` : Table de référence des 20 stations inventoriées (coordonnées, bassin, pas de temps, période, couverture 2016–2020, flags de qualité).
  - `obs_daily.csv` : Séries temporelles de débit journalier ($m^3/s$).
  - `obs_monthly.csv` : Séries temporelles de débit mensuel ($m^3/s$).
  - `README.txt` : Documentation complète sur l'origine des fichiers (ABH Sebou, ABH Bouregreg-Chaouia), les conversions géodésiques Lambert $\rightarrow$ WGS84 et les audits de qualité.
- **Répertoire des données brutes sources :**  
  `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/streamflow_data/`

---

### 2.2 Bassin du Bouregreg
- **Pas de temps natif :** **Journalier** (`obs_daily.csv`).
- **Période historique couverte :** De 1971/1975 jusqu'au 31 août 2023.
- **Couverture sur l'année pilote 2019 :** **365 jours complets (100 %)** pour toutes les stations instrumentées.
- **Détail des 6 stations :**
  1. `b_l_chafia` (**Lalla Chafia**, oued Bouregreg) :
     - Coordonnées : Lon $-6{,}3800^\circ$, Lat $+33{,}7000^\circ$ (source : lu_wgs84).
     - Période : 1971-08-16 au 2023-08-31 (19 007 valeurs, 99{,}95 % sur 2016–2020, 1 jour manquant le 2016-10-06 ; 365 j en 2019).
     - Statut HyMAP : Appariée et prioritaire.
  2. `b_s_jabeur` (**Sidi Jabeur**, oued Grou) :
     - Coordonnées : Lon $-6{,}4300^\circ$, Lat $+33{,}5800^\circ$ (source : lu_wgs84).
     - Période : 1971-11-07 au 2023-08-31 (18 925 valeurs, 100 % sur 2016–2020 ; 365 j en 2019).
     - Statut HyMAP : Appariée et prioritaire.
  3. `b_a_loudah` (**Ain Loudah**, oued Korifla) :
     - Coordonnées : Lon $-6{,}7600^\circ$, Lat $+33{,}5500^\circ$ (source : lu_wgs84).
     - Période : 1971-07-22 au 2023-08-31 (19 013 valeurs, 100 % sur 2016–2020 ; 365 j en 2019).
     - Statut HyMAP : Appariée (attention : doublon avec S.M. Cherif).
  4. `b_s_cherif` (**S.M. Cherif**, oued Mechraa) :
     - Coordonnées : Lon $-6{,}6300^\circ$, Lat $+33{,}5500^\circ$ (source : lu_wgs84).
     - Période : 1971-07-22 au 2023-08-31 (19 013 valeurs, 100 % sur 2016–2020 ; 365 j en 2019).
     - Statut HyMAP : Appariée (attention : doublon avec Ain Loudah).
  5. `b_a_ziar` (**Ain Ziar**) :
     - Données : 1975-10-02 au 2023-08-31 (17 500 valeurs, 365 j en 2019).
     - Statut HyMAP : **Inexploitable pour le routage** (coordonnées manquantes dans tous les fichiers sources ABH).
  6. `b_r_fathia` (**Ras Fathia**) :
     - Données : 1975-09-01 au 2023-08-31 (17 531 valeurs, 365 j en 2019).
     - Statut HyMAP : **Inexploitable pour le routage** (coordonnées manquantes dans tous les fichiers sources ABH).

> [!WARNING]
> **Anomalie critique de duplication Ain Loudah / S.M. Cherif :**  
> Les séries journalières observées d'Ain Loudah et S.M. Cherif sont identiques à $96{,}3\,\%$ sur l'ensemble de l'historique et rigoureusement identiques ($100\,\%$) sur 2016–2020 (1 827 jours identiques). Ces deux enregistrements ne constituent pas deux contrôles indépendants. Le bassin du Bouregreg compte en réalité **3 stations physiques indépendantes exploitables pour HyMAP** (`b_l_chafia`, `b_s_jabeur`, et la série commune `b_a_loudah` / `b_s_cherif`).

---

### 2.3 Bassin du Sebou
- **Pas de temps natif :** **Mensuel** (`obs_monthly.csv`). **Aucune mesure journalière n'est disponible sur le cluster pour le Sebou.**
- **Période historique couverte :** De 1935 à août 2024 selon les stations.
- **Couverture sur l'année pilote 2019 :** Variable selon les stations (de 0 à 12 mois).
- **Détail des 11 stations :**
  1. `ain_timedrine` (**Ain Timedrine**, oued Sebou, IRE 581/22) :
     - Coordonnées : Lon $-4{,}54864^\circ$, Lat $+33{,}74755^\circ$, Aire : $4\,387\text{ km}^2$.
     - Période : 1935-02 à 2024-08 (830 mois, 100 % sur 2016–2020 ; **12 mois en 2019**).
     - Statut HyMAP : Appariée.
  2. `azzaba` (**Azzaba**, oued Sebou, IRE 583/22) :
     - Coordonnées : Lon $-4{,}64748^\circ$, Lat $+33{,}83502^\circ$, Aire : $4\,640\text{ km}^2$.
     - Période : 1957-09 à 2024-08 (804 mois, 100 % sur 2016–2020 ; **12 mois en 2019**).
     - Statut HyMAP : Appariée.
  3. `dar_el_hamra` (**Dar El Hamra**, oued Zloul, IRE 1000/23) :
     - Coordonnées : Lon $-4{,}41204^\circ$, Lat $+33{,}76728^\circ$, Aire : $670\text{ km}^2$.
     - Période : 1984-08 à 2024-08 (481 mois, 100 % sur 2016–2020 ; **12 mois en 2019**).
     - Statut HyMAP : Appariée.
  4. `pont_ait_aissa` (**Pont Ait Aissa**, oued Guigou, IRE 198/30) :
     - Coordonnées : Lon $-5{,}06952^\circ$, Lat $+33{,}16697^\circ$, Aire : $205\text{ km}^2$.
     - Période : 1999-09 à 2024-08 (300 mois, 100 % sur 2016–2020 ; **12 mois en 2019**).
     - Statut HyMAP : Appariée.
  5. `ait_khabbach` (**Ait Khabbach**, oued Guigou, IRE 585/22) :
     - Coordonnées : Lon $-4{,}78597^\circ$, Lat $+33{,}43057^\circ$, Aire : $1\,240\text{ km}^2$.
     - Période : 1970-09 à 2024-08 (591 mois, 73{,}3 % sur 2016–2020 ; **8 mois en 2019 : janvier à août**, lacune sep–déc).
     - Statut HyMAP : Appariée.
  6. `pont_de_mdez` (**Pont de Mdez**, oued Sebou, IRE 582/22) :
     - Coordonnées : Lon $-4{,}52235^\circ$, Lat $+33{,}67483^\circ$, Aire : $3\,435\text{ km}^2$.
     - Période : 1951-09 à 2019-08 (816 mois, 73{,}3 % sur 2016–2020 ; **8 mois en 2019 : janvier à août**, arrêt définitif en août 2019).
     - Statut HyMAP : Appariée.
  7. `el_mers` (**El Mers**, oued Masser, IRE 541/23) :
     - Coordonnées : Lon $-4{,}39829^\circ$, Lat $+33{,}46408^\circ$, Aire : $1\,004\text{ km}^2$.
     - Période : 1981-09 à 2017-08 (456 mois ; **0 mois en 2019**, arrêt en août 2017).
     - Statut HyMAP : Appariée historique (inutilisable en 2019).
  8. **Stations historiques sans données récentes (arrêtées avant 2016) :**
     - `ain_aicha` (Ouergha, IRE 1217/9, $2\,420\text{ km}^2$) : arrêt en août 2015.
     - `ain_al_ouali` (Lyhoudi, IRE 2210/15, $408\text{ km}^2$) : arrêt en août 2015.
     - `azib_soltane` (Sebou, IRE 1540/15, $16\,150\text{ km}^2$) : arrêt en août 2015.
     - `allal_tazi` (Sidi Allal Tazi, aval Sebou, IRE 1355/8, $26\,600\text{ km}^2$) : arrêt en septembre 1991.

---

### 2.4 Synthèse des Stations Disponibles pour l'Évaluation 2019

| Bassin | Identifiant | Nom Station | Oued | Pas de Temps | Couverture 2019 | Échantillon 2019 | Statut et Remarque |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| **Bouregreg** | `b_l_chafia` | Lalla Chafia | Bouregreg | Journalier | 100 % | 365 jours | **Valide, indépendante, prioritaire** |
| **Bouregreg** | `b_s_jabeur` | Sidi Jabeur | Grou | Journalier | 100 % | 365 jours | **Valide, indépendante, prioritaire** |
| **Bouregreg** | `b_a_loudah` | Ain Loudah | Korifla | Journalier | 100 % | 365 jours | Valide (doublon strict S.M. Cherif) |
| **Bouregreg** | `b_s_cherif` | S.M. Cherif | Mechraa | Journalier | 100 % | 365 jours | Valide (doublon strict Ain Loudah) |
| **Bouregreg** | `b_a_ziar` | Ain Ziar | - | Journalier | 100 % | 365 jours | Inexploitable (coordonnées absentes) |
| **Bouregreg** | `b_r_fathia` | Ras Fathia | - | Journalier | 100 % | 365 jours | Inexploitable (coordonnées absentes) |
| **Sebou** | `ain_timedrine` | Ain Timedrine | Sebou | Mensuel | 100 % | 12 mois | **Valide (pas mensuel exclusivement)** |
| **Sebou** | `azzaba` | Azzaba | Sebou | Mensuel | 100 % | 12 mois | **Valide (pas mensuel exclusivement)** |
| **Sebou** | `dar_el_hamra` | Dar El Hamra | Zloul | Mensuel | 100 % | 12 mois | **Valide (pas mensuel exclusivement)** |
| **Sebou** | `pont_ait_aissa`| Pont Ait Aissa | Guigou | Mensuel | 100 % | 12 mois | **Valide (pas mensuel exclusivement)** |
| **Sebou** | `ait_khabbach` | Ait Khabbach | Guigou | Mensuel | 66{,}7 % | 8 mois | Partielle (janvier à août 2019) |
| **Sebou** | `pont_de_mdez` | Pont de Mdez | Sebou | Mensuel | 66{,}7 % | 8 mois | Partielle (s'arrête en août 2019) |
| **Sebou** | `el_mers` | El Mers | Masser | Mensuel | 0 % | 0 mois | Inexploitable en 2019 (arrêt 2017) |

---

## 3. Configuration et Scripts du Routage HyMAP2 de P1

### 3.1 Mode de Fonctionnement : Strictement Offline (Découplé)
Dans la chaîne de modélisation de la Phase 1, **le routage hydrodynamique HyMAP2 n'est pas couplé de façon interactive dans l'exécutable Noah-MP.**
Il fonctionne selon une chaîne séquentielle en 3 étapes distinctes :
1. **Simulation Noah-MP 4.0.1 (Surface) :**
   Noah-MP s'exécute au pas de temps de 15 minutes sur la grille régulière nodale $181 \times 121$ ($0{,}05^\circ$) et exporte quotidiennement dans ses fichiers NetCDF d'historique (`SURFACEMODEL/YYYYMM/LIS_HIST_YYYYMMDD0000.d01.nc`) les moyennes sur 24 heures de deux flux clés :
   - `Qs_tavg` : ruissellement de surface rapide ($kg \cdot m^{-2} \cdot s^{-1}$) ;
   - `Qsb_tavg` : drainage sous-racinaire profond de base ($kg \cdot m^{-2} \cdot s^{-1}$).
2. **Remapping Spatial des Ruissellements :**
   Un script Python reprojette au pas journalier les flux `Qs_tavg` et `Qsb_tavg` depuis la grille nodale $181 \times 121$ vers la grille centrée HyMAP $180 \times 120$ via un stencil arithmétique 4 points.
3. **Exécution du Routage HyMAP2 via LIS :**
   Le binaire LIS compilé avec le modèle de routage HyMAP2 est exécuté de manière autonome (*offline*), forcé par les fichiers de ruissellement remappés. HyMAP propage les hydrogrammes en lit mineur et en nappe et écrit les fichiers d'historique de routage (`ROUTING/YYYYMM/LIS_HIST_YYYYMMDD0000.d01.nc`) contenant `Streamflow_tavg` ($m^3/s$).

---

### 3.2 Emplacement des Fichiers et Scripts sur le Cluster
- **Répertoire de travail HyMAP de P1 :**  
  `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/routing/`
- **Documentation canonique P1 :**
  - `README_HYMAP.md` (synthèse complète des configurations, spin-up et correctifs) ;
  - `docs/HYMAP2_ROUTING_GUIDE.md` (guide technique détaillé de 43 pages) ;
  - `docs/POSTPROCESSING.md` (protocole de post-traitement et reproduction).
- **Fichiers de paramètres physiographiques (Augusto Getirana, NASA GSFC) :**
  - Archive source MERIT-Hydro à $0{,}05^\circ$ : `data/Hymap_params_from_Augusto/lis_param.tar.gz`.
  - Grille assemblée par LDT : `ldt_hymap/lis_input_corrigee.d01.nc`.
  - **Correctif physico-numérique critique de nappe :** application obligatoire du diviseur 86 400 sur la constante de vidange `tbasflw` dans `ldt_hymap/ldt.config` (conversion des secondes fournies en jours attendus par LISF, évitant l'accumulation continue et irréaliste de stock souterrain).
- **Scripts de remapping du ruissellement :**
  - OPL : `scripts/remap_runoff_opl.py`
  - DA-NoCDF : `scripts/remap_runoff_da_nocdf.py`
  - DA-CDF : `scripts/remap_runoff_da_cdf.py`
  - Emplacement des ruissellements reprojetés : `runoff_remap_180x120/<EXP>/SURFACEMODEL/` (1 828 fichiers journaliers NetCDF par expérience).
- **Configurations LIS pour HyMAP2 (`configs/`) :**
  - Spin-up Cycle 1 OPL corrigée (2016) : `lis_hymap2_opl_spinup_c1_corrigee.config`
  - Spin-up Cycle 2 OPL corrigée (2016) : `lis_hymap2_opl_spinup_c2_corrigee.config`
  - Production 5 ans OPL corrigée (2016–2020) : `lis_hymap2_opl_2016_2020_corrigee.config`
  - Chaîne DA-NoCDF : `lis_hymap2_da_nocdf_spinup_c1.config`, `c2.config`, `2016_2020.config`
  - Chaîne DA-CDF : `lis_hymap2_da_cdf_spinup_c1.config`, `c2.config`, `2016_2020.config`
- **Scripts SLURM (`jobs/`) :**
  - `run_hymap2_opl_2016_2020_corrigee.sh`, `run_hymap2_da_nocdf_2016_2020.sh`, `run_hymap2_da_cdf_2016_2020.sh` (1 nœud, 16 cœurs MPI, ~12 min pour 1 an de spin-up, ~45 min pour 5 ans de production).
- **Extraction des débits et évaluation hydrologique :**
  - Table d'appariement des mailles HyMAP : `stations/hymap_stations_pairing.csv` (12 stations).
  - Script d'extraction : `scripts/extract_streamflow.py` (extrait la variable `Streamflow_tavg` en $m^3/s$).
  - Script d'évaluation : `scripts/evaluate_sebou_bouregreg.py` (calcule KGE, NSE, $r$, PBIAS contre `obs_daily.csv` et `obs_monthly.csv`).
  - Caches journaliers de référence :
    - OPL corrigée : `output/OPL_2016_2020_corrigee/t7_sim_daily_12stations_corrigee.csv`
    - DA-NoCDF : `output/DA_nocdf_2016_2020/t8_sim_daily_12stations_da_nocdf.csv`
    - DA-CDF : `output/DA_cdf_2016_2020/t9_sim_daily_12stations_da_cdf.csv`
- **Restarts de nappe et rivière disponibles au 2019-01-01 :**
  - Des fichiers de restart `LIS_RST_HYMAP2_router_201812312345.d01.nc` existent dans les répertoires `output/<EXP>/ROUTING/201812/` des simulations 5 ans de P1 (OPL, DA-CDF, DA-NoCDF).

---

## 4. Produits d'Évapotranspiration (ET) Disponibles

| Produit | Variable(s) | Résolution Spatiale | Pas de Temps | Période sur le Cluster | Couverture 2019 | Chemin d'Accès Cluster | Statut / Utilisation P1 |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- |
| **GLEAM v4.3a** (Global Land Evaporation Amsterdam Model) | $E$ (ET réelle), $E_p$ (ET potentielle) | $0{,}25^\circ$ ($pprox 25\text{ km}$) | Journalier | 2015–2020 complet | 100 % (365 j) | `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/evapotranspiration/GLEAM/raw/` | **Référence principale P1** ; harmonisée et masquée sur 4 110 mailles continentales dans `smap-xai-audit/data/processed/value/gleam/reference/` |
| **MODIS MOD16A2GF** (Collection 6.1) | `ET_500m`, `PET_500m`, `ET_QC_500m` | 500 m | Composite 8 jours | 2015–2020 complet (46 fichiers/an) | 100 % (46 dates) | `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/evapotranspiration/MOD16/raw/` | **Complet** ; tuile sinusoïdale `h17v05` couvrant l'intégralité du Nord Maroc |
| **WaPOR v3 AETI** (FAO Water Productivity) | `AETI` (ET réelle annuelle / mensuelle) | ~250 m | Annuel (fichiers dispos) | 2018, 2019, 2020 | 100 % (annuel) | `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/evapotranspiration/WaPOR/WaPOR_v3_AETI/` | Rasters GeoTIFF disponibles ; niveau continental africain |
| **GLDAS CLSM 1.0** (NASA GSFC) | `Evap_tavg` ($kg \cdot m^{-2} \cdot s^{-1}$) | $1{,}0^\circ$ | 3-horaire | 2015–2020 complet | 100 % | `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/lsm/GLDAS/GLDAS_CLSM10_3H/` | Modèle d'intercomparaison continentale globale |
| **GLDAS NOAH 0.25** (NASA GSFC) | `Evap_tavg` ($kg \cdot m^{-2} \cdot s^{-1}$) | $0{,}25^\circ$ | 3-horaire | 2015 à juin 2018 | **0 % en 2019** | `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/lsm/GLDAS/GLDAS_NOAH025_3H/` | Partiel (interrompu en juin 2018) |

---

## 5. Mesures In Situ d'Humidité du Sol et Produits de Substitution

### 5.1 Constat Formel sur les Mesures In Situ
- **Aucune mesure in situ d'humidité du sol (stations agrométéorologiques, sondes TDR, sondes capacitives, profilomètres, stations de validation ISMN ou SMAP Cal/Val) n'est disponible sur le cluster pour le Maroc.**
- Le répertoire réservé `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/soil_moisture/in_situ/` est strictement vide.
- Les fichiers du répertoire `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/data/insitu_data/` (`Chefchaouen.txt`, `Agadir-massira.txt`, `ABH_stations_plusziz_and_abhbc1.csv`) contiennent des relevés pluviométriques et thermométriques de surface des stations météorologiques de la Direction Générale de la Météorologie (DGM) et des Agences de Bassin, sans aucune observation d'humidité du sol.
- L'audit méthodologique de P1 (`smap-xai-audit/docs/nocdf_native_rf_results.md`) a formellement acté cette absence d'observations in situ pour l'humidité du sol à l'échelle du domaine d'étude.

---

### 5.2 Produits Satellitaires d'Humidité du Sol Indépendants Disponibles
Pour compenser l'absence de mesures in situ, deux jeux de données satellitaires micro-ondes indépendants de SMAP sont présents sur le cluster :

1. **Copernicus Global Land Service - ASCAT Surface Soil Water Index (SWI) :**
   - **Emplacement :**  
     `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/soil_moisture/ASCAT/raw/`
   - **Versions :** v3.1.1 et v4.0.1 (fichiers journaliers globaux NetCDF).
   - **Période :** 2015–2020 complet (couverture intégrale 2019).
   - **Résolution :** Grille globale $0{,}1^\circ$ ($pprox 10\text{ km}$), pas de temps journalier.
   - **Intérêt pour l'humidité racinaire :** Le produit fournit le Soil Water Index pour différents facteurs de filtrage exponentiel $T$ ($T=1, 5, 10, 15, 20, 40, 60, 100\text{ jours}$). Les valeurs $T \ge 20$ intègrent la mémoire temporelle d'infiltration et constituent un proxy reconnu de l'humidité de la zone racinaire (*root-zone soil moisture*, RZSM).
2. **ESA CCI Soil Moisture Combined (v08.1) :**
   - **Emplacement :**  
     `/home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/data/validation/soil_moisture/ESA_CCI/raw/`
   - **Période :** 2015–2020 complet (couverture intégrale 2019).
   - **Résolution :** Grille régulière $0{,}25^\circ$, pas de temps journalier.
   - **Contenu :** Produit multi-capteurs combiné actif/passif d'humidité superficielle ($m^3/m^3$).

---

## 6. Synthèse d'Adéquation Technique pour les Critères de la Décision 16

1. **Critère C1 ($|\text{incrément moyen couche 1}|$ sur le domaine, DA-V2 < DA-NoCDF) :**  
   - **Faisabilité : 100 %**. Calcul direct sur les sorties NetCDF journalières des trois chaînes 2019 (`output/<EXP>/SURFACEMODEL/` ou fichiers de diagnostic EnKF).
2. **Critère C2 ($|\text{eau cumulée ajoutée par l'assimilation en 2019}|$ en mm, moyenne domaine, DA-V2 < DA-NoCDF) :**  
   - **Faisabilité : 100 %**. Calcul direct sur les incréments EnKF intégrés verticalement sur les 4 couches de sol pondérées par leurs épaisseurs ($0{,}1\text{ m}$, $0{,}3\text{ m}$, $0{,}6\text{ m}$, $1{,}0\text{ m}$).
3. **Critère C3 (KGE des débits HyMAP aux stations du Sebou et du Bouregreg, DA-V2 $\ge$ KGE(DA-CDF-train) $- 0{,}05$) :**  
   - **Faisabilité partielle et contraintes méthodologiques majeures :**
     - **Bassin du Bouregreg :** Faisable au pas de temps **journalier** sur les 3 stations physiques indépendantes (`b_l_chafia`, `b_s_jabeur`, `b_a_loudah`), qui disposent de 365 observations en 2019.
     - **Bassin du Sebou :** **Non faisable au pas journalier** (aucune donnée journalière n'existe). Au pas **mensuel**, seules 4 stations disposent des 12 mois de 2019 (`ain_timedrine`, `azzaba`, `dar_el_hamra`, `pont_ait_aissa`), 2 stations s'arrêtent en août 2019 (8 mois), et `el_mers` s'est arrêtée en 2017 (0 mois). Calculer un KGE sur $N=12$ points mensuels présente une fragilité statistique extrême.
     - **Chaîne de calcul HyMAP :** HyMAP2 étant un routage *offline*, il nécessite l'exécution d'une chaîne post-traitement dédiée (remapping spatial des ruissellements $Q_s$ et $Q_{sb}$, puis exécution de LIS HyMAP2).
     - **Initialisation hydraulique au 2019-01-01 :** DA-V2 ne disposant pas d'une simulation 2016–2018 pour réaliser le spin-up HyMAP spécifique en 2 cycles, le routage 2019 doit obligatoirement être initialisé soit par un restart donneur de P1 (ex: `LIS_RST_HYMAP2_router_201812312345.d01.nc` issu de DA-CDF), soit en coldstart (ce qui vicie les 2 à 3 premiers mois de 2019).
