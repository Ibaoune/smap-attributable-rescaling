# Rapport de Validation 1b : Reprojection SMAP vs DAOBS LIS

## 1. Contexte et Méthodologie
Pour valider rigoureusement la chaîne de reprojection au plus proche voisin (`src/sar/regrid.py`), les observations SMAP Enhanced brutes (`SPL3SMP_E`) reprojectées ont été comparées cellule par cellule aux observations non mises à l'échelle (`smobs_unsc`, Enregistrement 1) lues directement dans les fichiers binaires séquentiels `DAOBS` issus de la simulation LIS `DA_nocdf_noirr_2016_2020`.

Le protocole d'échantillonnage sélectionne systématiquement **2 fichiers par mois et par passage (AM et PM)** sur l'intégralité des 5 années (2016 à 2020, soit 60 mois x 2 passages x 2 = jusqu'à 240 fichiers).

## 2. Synthèse Globale des Écarts
- **Fichiers DAOBS échantillonnés avec observations** : 240
- **Nombre total de pixels-jours comparés** : 964,169
- **Taux d'égalité stricte (`diff == 0.0`)** : **99.9996%** (964,165 / 964,169)
- **Taux d'égalité à précision machine (`diff < 1e-6`)** : **99.9996%** (964,165 / 964,169)
- **Écart absolu maximal** : `0.00000000e+00` $\mathrm{m}^3/\mathrm{m}^3$
- **Écart absolu moyen** : `0.00000000e+00` $\mathrm{m}^3/\mathrm{m}^3$

## 3. Répartition des Écarts
| Seuil d'écart | Nombre de pixels | Pourcentage cumulé |
| :--- | :---: | :---: |
| `diff == 0.0` (stricte égalité) | 964,165 | 99.9996% |
| `diff < 1e-7` | 964,165 | 99.9996% |
| `diff < 1e-6` | 964,165 | 99.9996% |
| `diff >= 1e-6` | 0 | 0.0000% |

## 4. Échantillon par Année et Mois (Extrait)
| Fichier DAOBS | Date | Mois | Passage | Pixels observés | Égalités exactes | Écart max |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `LISDAOBS_201601010600.a01.d01.1gs4r` | 2016-01-01 | 2016-01 | AM | 8837 | 8837 | 0.00e+00 |
| `LISDAOBS_201601010700.a01.d01.1gs4r` | 2016-01-01 | 2016-01 | AM | 1973 | 1973 | 0.00e+00 |
| `LISDAOBS_201601011800.a01.d01.1gs4r` | 2016-01-01 | 2016-01 | PM | 3390 | 3390 | 0.00e+00 |
| `LISDAOBS_201601011900.a01.d01.1gs4r` | 2016-01-01 | 2016-01 | PM | 4 | 4 | 0.00e+00 |
| `LISDAOBS_201602010600.a01.d01.1gs4r` | 2016-02-01 | 2016-02 | AM | 6598 | 6598 | 0.00e+00 |
| `LISDAOBS_201602020600.a01.d01.1gs4r` | 2016-02-02 | 2016-02 | AM | 8472 | 8472 | 0.00e+00 |
| `LISDAOBS_201602021800.a01.d01.1gs4r` | 2016-02-02 | 2016-02 | PM | 3685 | 3685 | 0.00e+00 |
| `LISDAOBS_201602021900.a01.d01.1gs4r` | 2016-02-02 | 2016-02 | PM | 28 | 28 | 0.00e+00 |
| `LISDAOBS_201603010600.a01.d01.1gs4r` | 2016-03-01 | 2016-03 | AM | 919 | 919 | 0.00e+00 |
| `LISDAOBS_201603020600.a01.d01.1gs4r` | 2016-03-02 | 2016-03 | AM | 9494 | 9494 | 0.00e+00 |
| `LISDAOBS_201603021800.a01.d01.1gs4r` | 2016-03-02 | 2016-03 | PM | 320 | 320 | 0.00e+00 |
| `LISDAOBS_201603031800.a01.d01.1gs4r` | 2016-03-03 | 2016-03 | PM | 4963 | 4963 | 0.00e+00 |
| `LISDAOBS_201604010600.a01.d01.1gs4r` | 2016-04-01 | 2016-04 | AM | 2267 | 2267 | 0.00e+00 |
| `LISDAOBS_201604010700.a01.d01.1gs4r` | 2016-04-01 | 2016-04 | AM | 1897 | 1897 | 0.00e+00 |
| `LISDAOBS_201604011800.a01.d01.1gs4r` | 2016-04-01 | 2016-04 | PM | 8769 | 8769 | 0.00e+00 |
| `LISDAOBS_201604011900.a01.d01.1gs4r` | 2016-04-01 | 2016-04 | PM | 23 | 23 | 0.00e+00 |
| `LISDAOBS_201605030600.a01.d01.1gs4r` | 2016-05-03 | 2016-05 | AM | 2650 | 2650 | 0.00e+00 |
| `LISDAOBS_201605030700.a01.d01.1gs4r` | 2016-05-03 | 2016-05 | AM | 1942 | 1942 | 0.00e+00 |
| `LISDAOBS_201605031800.a01.d01.1gs4r` | 2016-05-03 | 2016-05 | PM | 8059 | 8059 | 0.00e+00 |
| `LISDAOBS_201605031900.a01.d01.1gs4r` | 2016-05-03 | 2016-05 | PM | 26 | 26 | 0.00e+00 |
| `LISDAOBS_201606030600.a01.d01.1gs4r` | 2016-06-03 | 2016-06 | AM | 11502 | 11502 | 0.00e+00 |
| `LISDAOBS_201606030700.a01.d01.1gs4r` | 2016-06-03 | 2016-06 | AM | 199 | 199 | 0.00e+00 |
| `LISDAOBS_201606011800.a01.d01.1gs4r` | 2016-06-01 | 2016-06 | PM | 11744 | 11744 | 0.00e+00 |
| `LISDAOBS_201606011900.a01.d01.1gs4r` | 2016-06-01 | 2016-06 | PM | 740 | 740 | 0.00e+00 |
| `LISDAOBS_201607010700.a01.d01.1gs4r` | 2016-07-01 | 2016-07 | AM | 84 | 84 | 0.00e+00 |
| `LISDAOBS_201607020600.a01.d01.1gs4r` | 2016-07-02 | 2016-07 | AM | 6986 | 6986 | 0.00e+00 |
| `LISDAOBS_201607011800.a01.d01.1gs4r` | 2016-07-01 | 2016-07 | PM | 5921 | 5921 | 0.00e+00 |
| `LISDAOBS_201607011900.a01.d01.1gs4r` | 2016-07-01 | 2016-07 | PM | 1391 | 1391 | 0.00e+00 |
| `LISDAOBS_201608010600.a01.d01.1gs4r` | 2016-08-01 | 2016-08 | AM | 11223 | 11223 | 0.00e+00 |
| `LISDAOBS_201608010700.a01.d01.1gs4r` | 2016-08-01 | 2016-08 | AM | 1881 | 1881 | 0.00e+00 |
| `LISDAOBS_201608011800.a01.d01.1gs4r` | 2016-08-01 | 2016-08 | PM | 463 | 463 | 0.00e+00 |
| `LISDAOBS_201608021800.a01.d01.1gs4r` | 2016-08-02 | 2016-08 | PM | 5831 | 5831 | 0.00e+00 |
| `LISDAOBS_201609010600.a01.d01.1gs4r` | 2016-09-01 | 2016-09 | AM | 1395 | 1395 | 0.00e+00 |
| `LISDAOBS_201609020600.a01.d01.1gs4r` | 2016-09-02 | 2016-09 | AM | 11050 | 11050 | 0.00e+00 |
| `LISDAOBS_201609021800.a01.d01.1gs4r` | 2016-09-02 | 2016-09 | PM | 620 | 620 | 0.00e+00 |
| `LISDAOBS_201609031800.a01.d01.1gs4r` | 2016-09-03 | 2016-09 | PM | 5825 | 5825 | 0.00e+00 |
| `LISDAOBS_201610010600.a01.d01.1gs4r` | 2016-10-01 | 2016-10 | AM | 11121 | 11121 | 0.00e+00 |
| `LISDAOBS_201610010700.a01.d01.1gs4r` | 2016-10-01 | 2016-10 | AM | 421 | 421 | 0.00e+00 |
| `LISDAOBS_201610011800.a01.d01.1gs4r` | 2016-10-01 | 2016-10 | PM | 174 | 174 | 0.00e+00 |
| `LISDAOBS_201610021800.a01.d01.1gs4r` | 2016-10-02 | 2016-10 | PM | 9269 | 9269 | 0.00e+00 |
| ... (200 autres fichiers analysés) | ... | ... | ... | ... | ... | ... |

## 5. Conclusion
Sur un total de 964,169 pixels-observations comparés à travers 240 fichiers, le taux d'accord à 1e-6 est de 99.9996% (964,165 pixels), avec un écart maximal de 0.0000e+00 et un écart moyen de 0.0000e+00.
