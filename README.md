# SMAP Attributable Rescaling (SAR)

## Objectif Scientifique
Ce dépôt développe et évalue une méthode de **correction de biais explicable et attribuable**
pour l'assimilation de l'humidité du sol satellitaire (SMAP) dans le modèle de surface continentale
Noah-MP (LIS-Hydro), en alternative au CDF matching empirique classique.

## Cible Principale
La variable cible est l'écart signé entre l'observation SMAP brute et l'humidité de surface de l'Open Loop :
$$d = \text{SMAP}_{\text{brut}} - \text{SSM}_{\text{OL}} \quad [\mathrm{m}^3/\mathrm{m}^3]$$
- **Couche 1** (0–10 cm de sol).
- **Passages AM et PM traités séparément** (matin ~06:00 locale, soir ~18:00 locale).
- **Moyenne journalière de l'OL** (`SoilMoist_tavg` 24h) retenue conformément à la décision actée.

## Domaine et Période
- **Domaine** : NorthMor (121 × 181 mailles à 0,05° de résolution spatiale, 16 030 pixels terrestres actifs).
- **Période** : 2016–2020 (5 années complètes).
- **Protocole d'évaluation** : Entraînement 2016–2019, année 2020 scellée (test indépendant), 5 plis spatiaux réguliers (blocs de 1° × 1°).

## Statut du Projet
- **Phase 0** : Initialisation du squelette, des configurations, du contrat de données et des modules réutilisables terminée.
- **Phase 1** : En cours de validation au Point de Contrôle 1.

## Structure du Dépôt
- `configs/` : Définitions déclaratives des chemins, cibles, prédicteurs et splits spatio-temporels.
- `docs/` : Périmètre (`scope.md`), contrat de données (`data_contract.md`), registre des décisions (`decisions.md`), et manifeste de provenance (`provenance_manifest.csv`).
- `src/sar/` : Bibliothèque modulaire réutilisable (`io`, `smap`, `regrid`, `qc`, `target`, `predictors`).
- `scripts/` : Scripts de traitement numérotés et reproductibles.
- `jobs/slurm/` : Scripts de soumission Slurm pour calculs lourds sur HPC.
- `data/` : Répertoires de données (`raw_copied/`, `interim/`, `processed/`, ignorés par git).
- `figures/`, `reports/` : Résultats graphiques et rapports de synthèse.
