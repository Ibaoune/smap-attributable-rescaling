# Slurm Job Submission Scripts

Ce répertoire héberge les scripts batch Slurm pour l'exécution des calculs lourds sur le cluster HPC :
- `01_copy_existing.sbatch`
- `02_regrid_smap.sbatch`
- `03_validate_daobs.sbatch`
- `04_derive_predictors.sbatch`
- `05_assemble_tables.sbatch`

Directives Slurm types :
- Partition : `standard` / `defq` (à adapter selon les files d'attente du cluster)
- CPUs : 8 à 16 threads
- Mémoire : 32 à 64 GB
- Environnement conda : `smap-xai-audit` ou `smap-attributable-rescaling`
