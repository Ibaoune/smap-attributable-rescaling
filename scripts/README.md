# Pipeline Scripts for SMAP Attributable Rescaling (SAR)

Les scripts sont numérotés séquentiellement pour garantir une reproductibilité stricte :

- `01_copy_existing_data.py` : Copie des données préexistantes depuis `smap-xai-audit` (OL, forçages, statiques, diagnostics EnKF) et inscription automatique dans `docs/provenance_manifest.csv`.
- `02_regrid_smap_raw.py` : Extraction et reprojection des fichiers SMAP Enhanced bruts (`SPL3SMP_E`) au plus proche voisin LIS sur la grille NorthMor (AM et PM).
- `03_validate_smap_vs_daobs.py` : Comparaison et validation de la reprojection SMAP face aux observations lues dans les fichiers DAOBS du run `DA_nocdf_noirr_2016_2020`.
- `04_derive_new_predictors.py` : Dérivation des prédicteurs physiques et dynamiques (MSWEP, IRRIGFRAC, SLOPETYPE, pente continue, cumuls de pluie, jours secs, indices de texture SOILPARM.TBL, LAI, saison).
- `05_assemble_modeling_tables.py` : Assemblage des tables de modélisation annuelles Parquet (2016 à 2020) intégrant la cible signée, les prédicteurs, les masques QC, et les plis spatio-temporels.
