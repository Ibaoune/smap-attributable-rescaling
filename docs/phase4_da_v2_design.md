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
