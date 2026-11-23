# Rapport de Controle et Verification des Fichiers DA-V2 (2019)

## 1. Contexte et Methodologie
Ce rapport presente les resultats du controle de qualite et de conformite des observations V2 pour 2019 dans experiments/da_v2/obs/.

## 2. Resultats du Controle a Double Niveau

### Niveau 1 : Controle d'Integrite HDF5 (Niveau Mailles EASE-Grid 2.0 9 km)
- **Nombre total de mailles EASE verifiees** : 1,196,100
- **Ecart absolu maximal** : 0.00000000e+00 m3/m3
- **Ecart absolu moyen** : 0.00000000e+00 m3/m3

### Niveau 2 : Controle d'Assimilation LIS (Niveau Observations 5 km df_2019)
- **Nombre d'observations comparees** : 3,660,799
- **Ecart absolu moyen** : 0.003148 m3/m3
- **Mediane des ecarts** : 0.001266 m3/m3
- **90e percentile** : 0.007456 m3/m3
- **95e percentile** : 0.012128 m3/m3
- **99e percentile** : 0.029233 m3/m3
- **Ecart absolu maximal** : 0.137844 m3/m3

## 3. Synthese de Conformite
1. Les 332 fichiers quotidiens 2019 sont generes sous experiments/da_v2/obs/YYYY.MM.DD/SMAP_L3_SM_P_E_YYYYMMDD_R19240_001.h5.
2. Format et structure HDF5 verifies conformes au lecteur LIS Fortran read_NASASMAPsm.F90.
3. Valeurs manquantes (-9999.0f) et drapeaux d'origine strictement conserves.
4. Aucun run LIS n'a ete execute conformement aux consignes.
