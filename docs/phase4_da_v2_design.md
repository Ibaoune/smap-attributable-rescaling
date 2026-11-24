# Spécification et Conception du Pilote d'Assimilation DA-V2 (2019)

## 1. Contexte et Objectifs
Ce document détaille la conception technique du pilote d'assimilation de données **DA-V2** sur l'année 2019. L'objectif est de fournir au système **LIS** (*Land Information System Framework* v7.5) des observations d'humidité du sol SMAP préalablement corrigées selon la méthode candidate **V2** (retrait des biais attribuables aux composantes amovibles et de contexte de représentation), avec la remise à l'échelle interne de LIS désactivée (`"none"`), **sans modifier aucune ligne de code source de LIS**.

Conformément au protocole, aucun run LIS n'est exécuté durant cette étape de préparation.

---

## 2. Analyse de la Configuration DA de Phase 1 (P1)

### 2.1. Lecteur LIS et Chaîne d'Assimilation SMAP
Dans la configuration Phase 1 (P1) de la matrice expérimentale NorthMor (`DA_cdf_noirr_2016_2020` et `DA_nocdf_noirr_2016_2020`), LIS lit et assimile SMAP via les modules suivants :
- **Plugin DA** : `src/lisf/lis/plugins/LIS_lsmda_pluginMod.F90` (lignes 2800–2812). Le couplage pour `Noah-MP.4.0.1` et `SMAP(NASA) soil moisture` (`NASASMAPsmobsId`) enregistre les routines :
  - Initialisation : `NoahMP401_dasoilm_init`
  - Extraction de l'état modèle : `NoahMP401_getsoilm`
  - Projection de l'état : `NoahMP401_getsmpred`
  - Contrôle qualité modèle : `NoahMP401_qc_soilmobs`
  - Mise à jour de l'analyse EnKF : `NoahMP401_updatesoilm`
- **Module lecteur d'observation** : `src/lisf/lis/dataassim/obs/NASA_SMAPsm/NASASMAPsm_Mod.F90` et `read_NASASMAPsm.F90`.
- **Routine de lecture spécifique 9 km** : Pour la désignation `SPL3SMP_E`, LIS appelle `read_NASASMAP_E_data` (`read_NASASMAPsm.F90`, lignes 208–260 et 924–1186).

### 2.2. Format des Fichiers et Nomenclature
- **Format** : Fichiers binaires hiérarchiques HDF5 (`.h5`).
- **Grille native** : Grille globale cylindrique EASE-Grid 2.0 9 km (M09), dimensions fixées dans `NASASMAPsm_Mod.F90` (lignes 372–373) :
  $$\text{nc} = 3856 \quad \text{colonnes}, \quad \text{nr} = 1624 \quad \text{lignes}$$
  (En convention C / Python HDF5, le dataset a pour forme `(1624, 3856)`).
- **Arborescence et nommage** : Définis dans `read_NASASMAPsm.F90` (lignes 220–224) :
  ```fortran
  list_files = trim(smobsdir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
               trim(dd)//'/SMAP_L3_SM_P_E_' &
               //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
               trim(CRID)//'*.h5'
  ```
  LIS recherche chaque jour les fichiers correspondant au motif :
  `${smobsdir}/YYYY.MM.DD/SMAP_L3_SM_P_E_YYYYMMDD_${CRID}*.h5`
  où `CRID` est le Composite Release ID configuré dans `lis.config` (dans P1 : `"R19"`).

### 2.3. Groupes et Variables Lus par LIS
Dans `read_NASASMAP_E_data` (`read_NASASMAPsm.F90`, lignes 960–970 et 1024–1140), LIS extrait exclusivement les datasets suivants selon le passage orbital :
- **Passage descendant AM ('D')** (acquisition vers ~06:00 UTC) :
  - Groupe : `Soil_Moisture_Retrieval_Data_AM`
  - Humidité du sol : `soil_moisture` (float32, valeur manquante `-9999.0f`, unité $\mathrm{cm}^3/\mathrm{cm}^3 = \mathrm{m}^3/\mathrm{m}^3$)
  - Drapeau qualité : `retrieval_qual_flag` (uint16)
  - Eau dans la végétation : `vegetation_water_content` (float32, unité $\mathrm{kg}/\mathrm{m}^2$)
- **Passage ascendant PM ('A')** (acquisition vers ~18:00 UTC) :
  - Groupe : `Soil_Moisture_Retrieval_Data_PM`
  - Humidité du sol : `soil_moisture_pm` (float32, valeur manquante `-9999.0f`)
  - Drapeau qualité : `retrieval_qual_flag_pm` (uint16)
  - Eau dans la végétation : `vegetation_water_content_pm` (float32)

### 2.4. Options de Remise à l'Échelle (Scaling Strategy)
La stratégie de rescaling est pilotée par le paramètre `Data assimilation scaling strategy:` lu dans `lis.config` :
- **Configuration CDF matching (P1 `DA_cdf_noirr_2016_2020`)** :
  - Fichier de configuration : `.../generated/configs/lis_DA_SMAP_cdf_noirr_2019-01.config`, ligne 77 :
    ```
    Data assimilation scaling strategy:                   "CDF matching"
    ```
  - Code source : `read_NASASMAPsm.F90`, lignes 540–567. LIS appelle `LIS_rescale_with_CDF_matching` en chargeant les CDF mensuelles observées et simulées spécifiées aux lignes 87–88 (`SMAP(NASA) model CDF file:` et `SMAP(NASA) observation CDF file:`).
- **Configuration sans CDF (P1 `DA_nocdf_noirr_2016_2020`)** :
  - Fichier de configuration : `.../generated/configs/lis_DA_SMAP_nocdf_noirr_2019-01.config`, ligne 77 :
    ```
    Data assimilation scaling strategy:                   "none"
    ```
  - Code source : `read_NASASMAPsm.F90`, lignes 539–580. Lorsque `dascaloption == "none"`, l'ensemble du bloc conditionnel de rescaling est contourné. LIS assigne directement la valeur de `sm_current` (issue de l'observation après interpolation) au vecteur d'assimilation `obsl` :
    ```fortran
    obsl = LIS_rc%udef
    do r = 1, LIS_rc%obs_lnr(k)
       do c = 1, LIS_rc%obs_lnc(k)
          if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
             obsl(LIS_obs_domain(n, k)%gindex(c, r)) = sm_current(c, r)
          endif
       enddo
    enddo
    ```

### 2.5. Règles de Contrôle Qualité (QC) Appliquées par LIS
LIS applique deux niveaux successifs de contrôle qualité :
1. **QC Interne du Lecteur SMAP** (`read_NASASMAPsm.F90`, lignes 1152–1163) :
   - Rejet si le contenu en eau de la végétation dépasse le seuil physique :
     $$\text{vegetation\_water\_content} > 5.0\,\mathrm{kg}/\mathrm{m}^2 \implies \text{sm\_data} = \text{LIS\_rc\%udef}$$
   - Rejet si le bit 0 du drapeau de qualité n'indique pas la recommandation officielle NASA (SMAP recommended quality) :
     $$\text{ibits}(\text{retrieval\_qual\_flag}, 0, 1) \ne 0 \implies \text{sm\_data} = -9999.0$$
   - Activé par la clé de configuration (ligne 85 de `lis.config`) : `SMAP(NASA) soil moisture apply SMAP QC flags: 1`.
2. **QC Dynamique du Modèle de Surface Noah-MP** (`noahmp401_qc_soilmobs.F90`, lignes 190–276) :
   - Précipitation active au moment du passage : $\text{rainf\_obs} > 3\times 10^{-6}\,\mathrm{kg}/\mathrm{m}^2/\mathrm{s}$ ($0{,}0108\,\mathrm{mm}/\mathrm{h}$) $\implies$ rejet.
   - Sol gelé : $|\theta_1 - \theta_{1,\text{liq}}| > 10^{-4}$ ou température de sol/surface $T_1 \le 273{,}15\,\mathrm{K}$ $\implies$ rejet.
   - Couverture neigeuse : équivalent en eau $\text{SWE} > 10^{-3}\,\mathrm{m}$ ou fraction de neige $\text{SCA} > 10^{-4}$ $\implies$ rejet.
   - Type de couverture terrestre : forêt dense ($\text{vegt} \le 4$) $\implies$ rejet.
   - Fraction de végétation verte : $\text{shdfac} > 0{,}90$ $\implies$ rejet.
   - Queues physiques de distribution : observation trop proche de la saturation ou du point de flétrissement :
     $$\theta_{\text{sat}} - \text{smobs} < 0{,}02 \quad \text{ou} \quad \text{smobs} - \theta_{\text{wlt}} < 0{,}02 \implies \text{rejet}.$$

---

## 3. Fournir à LIS des Observations Déjà Corrigées (DA-V2) sans Modifier le Code LIS

Pour injecter la correction $\mathrm{SMAP}^*_{\mathrm{V2}}$ dans LIS sans toucher à une seule ligne de son code Fortran, la stratégie repose sur le principe de l'observation pré-corrigée (*offline calibrated observations*) :

### 3.1. Modifications de la Configuration LIS (`lis.config`)
En s'appuyant sur le modèle Phase 1 `lis_da.config.template` et `lis_DA_SMAP_nocdf_noirr_2019-01.config` :
1. **Désactivation du rescaling LIS** (ligne 77) :
   ```
   Data assimilation scaling strategy:                   "none"
   ```
   *Effet* : LIS n'applique aucun CDF matching ni rescaling linéaire/anomalie ; il assimile les valeurs fournies telles quelles.
2. **Redirection du répertoire d'observations** (ligne 84) :
   ```
   SMAP(NASA) soil moisture data directory:  ./experiments/da_v2/obs
   ```
   *Effet* : LIS cherche ses fichiers quotidiens dans le répertoire des observations pré-corrigées V2.
3. **Maintien des paramètres de lecteur SMAP** (lignes 83, 85, 86) :
   ```
   SMAP(NASA) soil moisture data designation:        SPL3SMP_E
   SMAP(NASA) soil moisture apply SMAP QC flags:                 1
   SMAP(NASA) soil moisture Composite Release ID: "R19"
   ```
   *Effet* : LIS utilise son lecteur HDF5 standard pour fichiers 9 km EASE-Grid 2.0 et applique le QC NASA (VWC $\le 5$, bit 0 = 0).
4. **Gestion des fichiers CDF** (lignes 87–91) :
   Lorsque `Data assimilation scaling strategy: "none"` est activé, LIS n'ouvre ni ne lit les fichiers CDF spécifiés par `SMAP(NASA) model CDF file:` et `SMAP(NASA) observation CDF file:`. Ces lignes peuvent pointer vers des fichiers factices ou être conservées sans effet.
5. **Initialisation par Restart OL** (ligne 105) :
   ```
   Noah-MP.4.0.1 restart file: /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/OPL_noirr_2016_2020/restarts/LIS_RST_NOAHMP401_201901010000.d01.nc
   ```

### 3.2. Spécification des Fichiers dans `experiments/da_v2/obs/`
- **Structure des répertoires** : `experiments/da_v2/obs/YYYY.MM.DD/` pour chaque date observée de 2019.
- **Nom du fichier** : `SMAP_L3_SM_P_E_YYYYMMDD_R19240_001.h5`
- **Format HDF5** :
  - Dimensions globales identiques à SMAP : $1624 \times 3856$.
  - Groupes : `Soil_Moisture_Retrieval_Data_AM` et `Soil_Moisture_Retrieval_Data_PM`.
  - Dans `Soil_Moisture_Retrieval_Data_AM` :
    - `soil_moisture` : contient $\mathrm{SMAP}^*_{\mathrm{V2}}$ sur les cellules du domaine observées et valides, et `-9999.0f` partout ailleurs (valeur manquante préservée).
    - `retrieval_qual_flag` : conserve strictement les entiers 16 bits d'origine.
    - `vegetation_water_content` : conserve strictement les valeurs réelles d'origine.
  - Dans `Soil_Moisture_Retrieval_Data_PM` :
    - `soil_moisture_pm` : contient $\mathrm{SMAP}^*_{\mathrm{V2}}$ (passage PM) et `-9999.0f` ailleurs.
    - `retrieval_qual_flag_pm` : drapeaux d'origine préservés.
    - `vegetation_water_content_pm` : VWC d'origine préservée.

---

## 4. Vérification du Fichier Restart OL au 2019-01-01
Le restart de la simulation Open-Loop (`OPL_noirr_2016_2020`) a été formellement localisé et vérifié sur le système de fichiers Lustre :
- **Chemin absolu vérifié** :
  ```
  /home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/experiments/NorthMor/matrix_2016_2020/OPL_noirr_2016_2020/restarts/LIS_RST_NOAHMP401_201901010000.d01.nc
  ```
- **Attributs du fichier** :
  - Taille : `1,5 Mo` (`1 559 864 octets`)
  - Date de création : `18 juillet 2026, 16:25 UTC`
  - Format : NetCDF-3 / NetCDF-4 conforme LIS Noah-MP 4.0.1
  - Contenu : États complets des 4 couches de sol (`SMC`, `SH2O`, `STC`), humidité de canopée, profil de température et variables de surface au 2019-01-01 00:00 UTC.

Ce fichier constitue la condition initiale exacte et non perturbée pour démarrer l'assimilation DA-V2 au 1er janvier 2019.

---

## 5. Synthèse des Fichiers et Lignes de Configuration Concernés

| Composant | Fichier Source / Cible | Lignes | Rôle / Paramètre |
| :--- | :--- | :--- | :--- |
| **Config DA** | `lis_DA_SMAP_nocdf_noirr_2019-01.config` | L77 | `Data assimilation scaling strategy: "none"` |
| **Config DA** | `lis_DA_SMAP_nocdf_noirr_2019-01.config` | L83–86 | Désignation `SPL3SMP_E`, répertoire `experiments/da_v2/obs`, flags QC `1`, CRID `"R19"` |
| **Config DA** | `lis_DA_SMAP_nocdf_noirr_2019-01.config` | L105 | Restart OL : `/.../OPL_noirr_2016_2020/restarts/LIS_RST_NOAHMP401_201901010000.d01.nc` |
| **Template DA** | `templates/lis_da.config.template` | L77, L84, L105 | Variables de substitution pour la génération des scripts d'expérience |
| **Code LIS (Lecteur)** | `src/lisf/lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90` | L220–224 | Recherche des fichiers HDF5 par date et CRID |
| **Code LIS (Lecteur)** | `src/lisf/lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90` | L1024–1140 | Lecture des groupes AM/PM et datasets `soil_moisture`, `retrieval_qual_flag`, `VWC` |
| **Code LIS (Lecteur)** | `src/lisf/lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90` | L1152–1163 | QC SMAP interne (VWC $\le 5$, bit 0 = 0) |
| **Code LIS (Lecteur)** | `src/lisf/lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90` | L1174–1180 | Interpolation plus proche voisin vers le domaine modèle via `neighbor_interp` |
| **Code LIS (Scaling)** | `src/lisf/lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90` | L539–580 | Dérivation du scaling : contourné quand `dascaloption == "none"` |
| **Code LIS (QC LSM)** | `src/lisf/lis/surfacemodels/land/noahmp.4.0.1/da_soilm/noahmp401_qc_soilmobs.F90` | L190–276 | QC dynamique Noah-MP (pluie, gel, neige, forêt, saturation, flétrissement) |

---

## 6. Analyse des Lecteurs d'Observation LIS sur Grille Modèle (Piste 5 km en Lecture Seule)

Conformément à la demande d'investigation en lecture seule, l'ensemble des lecteurs d'humidité du sol disponibles dans `src/lisf/lis/dataassim/obs/` a été analysé pour déterminer si LIS propose un lecteur acceptant des observations déjà projetées sur la grille domaine LIS (5 km NorthMor) sans modification du code source LIS.

### 6.1. Lecteur `syntheticsm` (Synthetic Soil Moisture)
- **Module et routines** : `src/lisf/lis/dataassim/obs/syntheticsm/syntheticsmobs_module.F90` et `read_syntheticsmobs.F90`.
- **Nom LIS / Plugin ID** : `Synthetic soil moisture` (`LIS_synsmId`).
- **Format attendu** : Fichiers NetCDF (ou binaire direct non formaté) nommés selon le motif :
  `${smobsdir}/YYYYMM/SimObs_YYYYMMDDHHMN.nc`
  Variable requise : `SoilMoist_tavg` (NetCDF-3 ou NetCDF-4).
- **Grille** : Directement sur la grille domaine LIS `(obs_lnc, obs_lnr)` sans sous-routine d'interpolation : `smobs(c, r)` est injecté directement dans le vecteur local `obsl(LIS_obs_domain(n,k)%gindex(c,r))` (`read_syntheticsmobs.F90`, lignes 125–136).
- **Contrôle Qualité (QC)** :
  - Vérification des valeurs manquantes (`obsl /= -9999.0`).
  - Seuils physiques stricts de remise à l'échelle : `MIN_SM_VALUE = 0.0001`, `MAX_SM_VALUE = 0.45` (`read_syntheticsmobs.F90`, ligne 60).
  - Contrôle qualité LSM via `call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"//trim(LIS_synsmId)//char(0), n, k, OBS_state)`.
- **Compatibilité avec la configuration P1 (Noah-MP 4.0.1)** : **INCOMPATIBLE sans modifier le code LIS**.
  Dans `src/lisf/lis/plugins/LIS_lsmda_pluginMod.F90`, le couplage DA pour `LIS_synsmId` n'est enregistré que pour les modèles de surface suivants :
  - `Noah.3.3` (lignes 1101–1120)
  - `Noah.3.6` (lignes 1679–1698)
  - `NoahMP.3.6` (lignes 2087–2106)
  - `CLSM-F2.5` (lignes 3257–3276)
  - `RUC.3.7` (lignes 3572–3591)
  Dans le bloc dédié à **Noah-MP 4.0.1** (`LIS_noahmp401Id`, lignes 2750–3253), **aucun appel d'enregistrement `registerlsmdainit` n'existe pour `LIS_synsmId`**. Définir `Data assimilation observation type: "Synthetic soil moisture"` dans `lis.config` avec Noah-MP 4.0.1 déclenche une interruption fatale d'exécution dans le gestionnaire ESMF (`LIS_lsmda_init` : clé `"Noah-MP.4.0.1+Synthetic soil moisture"` non trouvée dans la table des plugins DA).

### 6.2. Lecteur `multisynsmobs` (Multiple Synthetic Soil Moisture)
- **Module et routines** : `src/lisf/lis/dataassim/obs/multisynsmobs/multisynsmobs_Mod.F90` et `read_multisynsmobs.F90`.
- **Nom LIS / Plugin ID** : `Multiple synthetic soil moisture` (`LIS_multisynsmId`).
- **Format attendu** : Fichier binaire brut direct séquentiel `${smobsdir}/YYYYMMDDHHMN.d01.gs4r` contenant 4 couches successives sur la grille (`sm1`, `sm2`, `sm3`, `sm4`).
- **Grille** : Grille domaine native LIS.
- **QC appliqué** : Test binaire `sm /= -9999.0`.
- **Compatibilité avec P1 (Noah-MP 4.0.1)** : **INCOMPATIBLE sans modifier le code LIS**. Non enregistré pour Noah-MP 4.0.1.

### 6.3. Lecteur `pildas` (PILDAS Synthetic Soil Moisture)
- **Module et routines** : `src/lisf/lis/dataassim/obs/pildas/pildassmobs_module.F90` et `read_pildassmobs.F90`.
- **Nom LIS / Plugin ID** : `PILDAS soil moisture` (`LIS_pildassmobsId`).
- **Format attendu** : Fichier NetCDF-4 `PILDAS-1_synobs_O002V001_YYYYMMDD_120000.nc4` avec le dataset `SFMC_OBS`.
- **Grille** : Projection matricielle PILDAS.
- **QC appliqué** : Borne physique `0.0 < SFMC_OBS <= 0.45`, rejet si $>0.50$, QC LSM.
- **Compatibilité avec P1 (Noah-MP 4.0.1)** : **INCOMPATIBLE sans modifier le code LIS**. Enregistré pour Noah-3.3 et NoahMP-3.6, mais absent pour Noah-MP 4.0.1.

### 6.4. Lecteur `SMAP_E_OPLsm` (SMAP Enhanced Open-Loop Synthetic)
- **Module et routines** : `src/lisf/lis/dataassim/obs/SMAP_E_OPLsm/SMAPEOPLsm_Mod.F90` et `read_SMAPEOPLsm.F90`.
- **Nom LIS / Plugin ID** : `SMAP E OPL soil moisture` (`LIS_SMAPEOPLsmobsId`).
- **Format attendu** : Fichiers NetCDF `${smobsdir}/ARFS_SM_*YYYYMMDDTHH*.nc`.
- **Grille** : Grille cylindrique équidistante (`MAP_PROJECTION = 'EQUIDISTANT CYLINDRICAL'`), interpolée vers la grille LIS via `neighbor_interp`. N'accepte pas directement la grille native LIS.
- **QC appliqué** : Alarme temporelle stricte, bornes `[0.0001, 0.45]`, QC modèle.
- **Compatibilité avec P1 (Noah-MP 4.0.1)** : Enregistré pour Noah-MP 4.0.1 (`LIS_lsmda_pluginMod.F90`, lignes 2854–2872), mais n'élimine pas l'étape d'interpolation et requiert un formatage externe ARFS spécifique non standard.

### 6.5. Synthèse et Justification de la Solution Retenue
| Lecteur LIS | Format de fichier | Grille d'entrée | Binding Noah-MP 4.0.1 | Utilisable sans modifier le code LIS ? |
| :--- | :--- | :--- | :--- | :--- |
| `NASA_SMAPsm` (**P1**) | HDF5 (`.h5`) | EASE-Grid 2.0 9 km | **Oui** (L2795) | **OUI (Recommandé, solution retenue)** |
| `syntheticsm` | NetCDF (`.nc`) | Grille LIS 5 km | **Non** (aucun binding) | **NON** (exige modification de `LIS_lsmda_pluginMod.F90`) |
| `multisynsmobs` | Binaire direct (`.gs4r`) | Grille LIS 5 km | **Non** | **NON** |
| `pildas` | NetCDF-4 (`.nc4`) | Grille LIS / PILDAS | **Non** | **NON** |
| `SMAP_E_OPLsm` | NetCDF (`.nc`) | Cylindrique équidistante | **Oui** (L2855) | **NON** (impose interpolation `neighbor_interp`) |
| `CDF_Transfer_NASA_SMAPsm` | HDF5 (`.h5`) | EASE-Grid 2.0 9 km | **Oui** (L2815) | **NON** (impose CDF transfer dynamique) |

**Conclusion** : Sous la contrainte absolue « sans modifier le code LIS », l'utilisation du lecteur natif **`NASA_SMAPsm`** (format HDF5 EASE-Grid 2.0 9 km) est l'unique solution fonctionnelle pour la configuration Noah-MP 4.0.1 de Phase 1. La préservation de l'intégrité de l'attribution V2 est garantie en masquant systématiquement à `-9999.0f` toute maille EASE du domaine NorthMor n'ayant pas reçu de prédiction V2, éliminant tout risque d'assimilation résiduelle de données SMAP brutes.


---

## 7. Initialisation des Runs DA en Phase 1 (P1) et Transposition aux Expériences 2019

### 7.1. Historique d'Initialisation en Phase 1 (P1)
En Phase 1 (matrice `matrix_2016_2020`), l'initialisation des runs DA s'est appuyée sur une chaîne hiérarchique en trois étapes :
1. **Spin-up Déterministe (`SPINUP_OPL`, 2005–2015)** :
   - Simulation déterministe mono-membre (`Number of ensembles per tile: 1`).
   - Le restart final au 2015-01-01 (`LIS_RST_NOAHMP401_201501010000.d01.nc`, 1,6 Mo) contient $N_{\text{tiles}} = 16\,030$ tiles.
2. **Expansion du Restart Mono-Membre vers Ensemble (20 membres)** :
   - **Absence d'option native dans LIS ou LDT** : Ni LIS ni LDT ne possèdent d'option intégrée (dans `lis.config` ou `ldt.config`) pour démultiplier un état restart Noah-MP NetCDF mono-membre vers $N$ membres.
   - **Script d'expansion Python** : L'auteur du projet a développé un script ad hoc (`NoahMP_Morocco/experiments/NorthMor/SPINUP/ENS_OPL_2015/scripts/expand_restart.py`), répliquant exactement 20 fois chaque variable le long de la dimension `ntiles` ($16\,030 \times 20 = 320\,600$ tiles) via `np.repeat(data, 20, axis=axis)`, produisant le fichier NetCDF 20 membres `LIS_RST_NOAHMP401_201501010000_ens20.d01.nc` (90 Mo).
3. **Génération de la Dispersion d'Ensemble (`ENS_OPL_2015`, janv–mars 2015)** :
   - Simulation LIS de 3 mois sans assimilation (`Number of data assimilation instances: 0`), 20 membres (`Number of ensembles per tile: 20`), initialisée avec `Start mode: restart` et `Perturbations start mode: "coldstart"`.
   - Les perturbations stochastiques GMAO sur le forçage atmosphérique et l'état de surface brisent la symétrie dès le premier pas de temps, amenant la dispersion d'ensemble à un niveau physique stable ($\\sigma \\approx 0{,}051\\,\\mathrm{m}^3/\\mathrm{m}^3$ au 2015-03-31).
4. **Spin-up d'Assimilation SMAP (`SPINUP_DA`, avr–déc 2015)** :
   - Initialisé au 2015-04-01 avec l'état de surface issu de `ENS_OPL_2015` au 2015-03-31 23:45 (`Start mode: restart`, fichier `restarts/surf/LIS_RST_NOAHMP401_201504010000.d01.nc`, 21 Mo).
   - Mode perturbations en `coldstart` pour le premier mois (`Perturbations start mode: "coldstart"`, `Perturbations restart filename: "none"`) suite à l'incompatibilité du fichier DAPERT issu d'ENS_OPL (Job 7260790 écourté par crash EOF car absence d'enregistrements d'observation dans DAPERT), puis en mode `restart` pour mai–décembre 2015.
5. **Démarrage de la Matrice 2016–2020 P1** :
   - Au 2016-01-01, les runs DA (`DA_nocdf_noirr_2016_2020` et `DA_cdf_noirr_2016_2020`) ont démarré directement depuis les états préconditionnés par 9 mois de DA SMAP issus de `SPINUP_DA` :
     - `Noah-MP.4.0.1 restart file: .../SPINUP_DA/restarts/surf/LIS_RST_NOAHMP401_201601010000.d01.nc` (L105)
     - `Perturbations restart filename: .../SPINUP_DA/restarts/pert/LIS_DAPERT_201601010000.d01.bin` (L173)
     - `Start mode: restart` (L27) et `Perturbations start mode: "restart"` (L171).

### 7.2. Paramétrage DA, Perturbations et Erreur d'Observation dans P1 (Citations de Config)
Dans les fichiers de configuration de production P1 (ex. `lis_DA_SMAP_nocdf_noirr_2019-01.config`) :
- **Taille d'ensemble** :
  ```
  L43:  Number of ensembles per tile:           20
  ```
- **Algorithme d'assimilation et variables d'état** :
  ```
  L67:  Number of data assimilation instances:                1
  L68:  Data assimilation algorithm:                          "EnKF"
  L72:  Data assimilation set:                                "SMAP(NASA) soil moisture"
  L75:  Data assimilation number of observation types:        1
  L76:  Data assimilation output ensemble spread:             1
  L80:  Number of state variables:                            4
  ```
- **Routage** :
  ```
  L66:  Routing model:                          "none"
  ```
  *(Le routage HyMAP n'est pas activé dans P1).*
- **Perturbations stochastiques (Schéma GMAO)** :
  ```
  L167: Forcing perturbation algorithm:           "GMAO scheme"
  L168: State perturbation algorithm:             "GMAO scheme"
  L169: Observation perturbation algorithm:       "GMAO scheme"
  L175: Forcing perturbation frequency:           "1hr"
  L177: Forcing perturbation attributes file:     "./data/pert_package/forcing_pertattribs.txt"
  L179: State perturbation frequency:              "3hr"
  L181: State perturbation attributes file:       "./data/pert_package/noahmp_sm_pertattribs.txt"
  L183: Observation perturbation frequency:         "6hr"
  L185: Observation perturbation attributes file:    "./data/pert_package/smap_pertattribs.txt"
  ```
- **Erreur d'observation SMAP** :
  Définie dans `./data/pert_package/smap_pertattribs.txt` :
  ```
  #perttype std  std_max  zeromean  tcorr  xcorr ycorr ccorr
  SMAP(NASA) soil moisture
    0   0.04     2.5      1  43200   0    0    1.0
  ```
  Bruit additif gaussien (`perttype: 0`), écart-type d'erreur $\sigma_o = 0{,}04\,\mathrm{m}^3/\mathrm{m}^3$ (`std: 0.04`), seuil de troncature à $2{,}5\,\sigma$ (`std_max: 2.5`), moyenne nulle garantie (`zeromean: 1`), temps de décorrélation temporelle $43\,200\,\text{s} = 12\,\text{h}$ (`tcorr: 43200`). L'écart-type d'erreur d'observation utilisé dans la mise à jour EnKF par LIS est exactement $0{,}04\,\mathrm{m}^3/\mathrm{m}^3$.

### 7.3. Protocole d'Initialisation pour les Expériences DA 2019
Pour garantir une comparaison strictement équitable et éliminer tout artefact transitoire au 2019-01-01, les trois expériences d'assimilation 2019 (`DA_nocdf_2019`, `DA_cdftrain_2019` et `DA_v2_2019`) suivront le protocole historique éprouvé lors de la Phase 1 (`ENS_OPL_2015` suivi de `SPINUP_DA` en avril 2015) :

1. **Spin-up d'ensemble Open-Loop préalable (`ENS_OL_2018Q4`)** :
   - Un run d'ensemble sans assimilation (`Number of data assimilation instances: 0`) est préparé sur les 3 mois précédant l'assimilation, du 2018-10-01 au 2019-01-01.
   - Il démarre du restart OL mono-membre au 2018-10-01 issu de P1 (`OPL_noirr_2016_2020/restarts/LIS_RST_NOAHMP401_201810010000.d01.nc`, 16 030 tiles), répliqué à 20 membres ($320\,600$ tiles) via `expand_restart.py` (`LIS_RST_NOAHMP401_201810010000_ens20.d01.nc`).
   - Les perturbations GMAO (forçage 1h, état 3h, observation 6h) sont activées en mode coldstart (`Perturbations start mode: "coldstart"`, `Perturbations restart filename: "none"`), selon la recette exacte d'`ENS_OPL_2015`.
   - Ce préconditionnement de 3 mois permet à l'ensemble d'acquérir un étalement dynamique et une structure de covariances d'erreur d'ébauche physiquement cohérents avant le début de l'assimilation.

2. **Démarrage des trois runs DA 2019 au 2019-01-01** :
   - Les trois runs 2019 démarreront du restart de surface 20 membres généré au 2019-01-01 par `ENS_OL_2018Q4` (`LIS_RST_NOAHMP401_201901010000.d01.nc`) :
     ```
     Start mode:                             restart
     Noah-MP.4.0.1 restart file:             experiments/da_v2/ENS_OL_2018Q4/output/2018Q4/SURFACEMODEL/201901/LIS_RST_NOAHMP401_201901010000.d01.nc
     ```
   - **Perturbations en coldstart** : Tout comme pour `SPINUP_DA` en avril 2015, les perturbations d'assimilation démarreront en mode coldstart :
     ```
     Perturbations start mode:               "coldstart"
     Perturbations restart filename:         "none"
     ```
     Ce choix évite tout problème de compatibilité de structure d'enregistrement dans les fichiers binaires `LIS_DAPERT` (qui diffèrent entre un run sans instances DA et un run avec 1 instance DA SMAP) tout en assurant que la dispersion d'ensemble des états de surface est déjà mature dès $t=0$.
   - Les trois configurations (`DA_nocdf_2019`, `DA_cdftrain_2019` et `DA_v2_2019`) partagent exactement le même état initial, garantissant que tout écart de performance provient exclusivement de la correction d'observation appliquée.

### 7.4. Définition Technique de la Métrique de Dispersion `ensspread` dans LIS

Dans le code source de LIS (version LISF 7.5) :
- Bien que la routine `writeEnsembleSpread` ([enkf_Mod.F90#L868](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/src/lisf/lis/dataassim/algorithm/enkf/enkf_Mod.F90#L868)) indique dans son en-tête `! This routine writes the ensemble spread (standard deviation) of the model state vector.`, le calcul effectif délégué à la sous-routine `LIS_writevar_spread` ([LIS_historyMod.F90#L7973-L7980](file:///home/mohammad.elaabaribao/lustre/empowermed-ahl6xm8o7mg/users/mohammad.elaabaribao/NoahMP_Morocco/src/lisf/lis/core/LIS_historyMod.F90#L7973-L7980)) ne calcule ni la variance ni l'écart-type.
- En effet, le calcul de la variance et de l'écart-type (`std_v = sqrt(var_v/npatch - mean_v*mean_v)`) a été mis sous directive conditionnelle inactive `#if 0` (lignes 7949–7970).
- Le bloc actif (lignes 7973–7980) calcule strictement l'**étendue d'ensemble (range)** :
  ```fortran
  min_v = 1000000
  max_v = -1000000 
  do t=1,npatch
     if(gtmp1(count1).lt.min_v) min_v = gtmp1(count1)
     if(gtmp1(count1).gt.max_v) max_v = gtmp1(count1)
     count1 = count1 + 1
  enddo                      
  gtmp(c,r) = max_v - min_v
  ```
- Pour chaque maille du domaine et chaque jour, la valeur archivée sous `ensspread_Soil Moisture Layer 1_01` correspond donc exactement à $\max(\text{smc}_1) - \min(\text{smc}_1)$ entre les 20 membres d'ensemble sur l'humidité volumique du sol de la couche 1 ($\text{m}^3/\text{m}^3$).

### 7.5. Incident de Concurrence `SMAP_filelist.dat` et Règle d'Isolation Stricte

#### 1. Mécanisme de l'incident
Lors de la lecture des observations SMAP, LIS délègue la recherche des fichiers HDF5 journaliers à une routine C POSIX (`create_filelist.c`) appelée depuis `read_NASASMAPsm.F90` :
```fortran
rc = create_filelist(trim(list_files)//char(0), "SMAP_filelist.dat"//char(0))
open (ftn, file="./SMAP_filelist.dat", action='read', status='old', iostat=ierr)
```
Le nom du fichier temporaire `"./SMAP_filelist.dat"` est écrit en dur relativement au répertoire d'exécution courant (`CWD`). Lorsque les trois expériences (`DA_nocdf_2019`, `DA_cdftrain_2019`, `DA_v2_2019`) s'exécutaient simultanément dans le répertoire racine `NoahMP_Morocco`, les processus concurrents écrasaient mutuellement ce fichier partagé à chaque cycle journalier (00:00).

La vérification exhaustive menée sur les fichiers `DAOBS` (Enregistrement 1, observations brutes non remises à l'échelle reprojectées à 5 km) et consignée dans `reports/check_obs_source_2019.csv` a révélé :
- `DA_nocdf_2019` : 8 cycles contaminés en janvier (lecture intempestive des fichiers V2 de `DA_v2`).
- `DA_cdftrain_2019` : 0 cycle contaminé en janvier, mais crash au 4 février suite à une troncature de chemin relative lors de l'accès concurrent.
- `DA_v2_2019` : **62 cycles contaminés** en janvier (lecture intempestive des fichiers SMAP bruts de `DA_nocdf`/`DA_cdftrain`).

#### 2. Règle de décision et remédiation appliquée
Conformément à la règle méthodologique fixée à l'avance :
1. **Annulation immédiate** : Les trois chaînes et le job d'orchestration ont été immédiatement annulés (`scancel`).
2. **Archivage** : L'ensemble des sorties corrompues (logs, outputs, restarts, configs générées) a été déplacé dans `experiments/da_v2/archive_filelist_bug/`.
3. **Isolation stricte des répertoires de travail** : Le générateur de scripts de soumission `chain_da.py` a été adapté afin que chaque expérience dispose de son propre répertoire d'exécution dédié :
   ```bash
   RUN_DIR="{base_dir}/run"
   mkdir -p "$RUN_DIR"
   for d in data src arch shared domain; do
       if [ -d "/home/.../NoahMP_Morocco/$d" ]; then
           ln -sfn "/home/.../NoahMP_Morocco/$d" "$RUN_DIR/$d"
       fi
   done
   cd "$RUN_DIR"
   ```
   Chaque expérience écrit et lit ainsi son propre fichier `./SMAP_filelist.dat` de façon strictement étanche.
4. **Relance propre au 2019-01-01** : Les trois chaînes ont été redémarrées simultanément au 2019-01-01 depuis le restart d'ensemble sain `ENS_OL_2018Q4` (`LIS_RST_NOAHMP401_201901010000.d01.nc`), perturbations en mode coldstart (`Perturbations restart filename: "none"`), sur la partition `compute` et le compte `empowermed-ahl6xm8o7mg-DEFAULT-CPU`.
