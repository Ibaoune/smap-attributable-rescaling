# Rapport de Validation 1b : Reprojection SMAP vs DAOBS LIS

## 1. Contexte et Méthodologie
Pour valider la chaîne de reprojection au plus proche voisin (`src/sar/regrid.py`), les observations SMAP Enhanced brutes (`SPL3SMP_E`) reprojectées ont été comparées cellule par cellule aux observations non mises à l'échelle (`smobs_unsc`) lues directement dans les fichiers binaires séquentiels `DAOBS` issus de la simulation LIS `DA_nocdf_noirr_2016_2020`.

## 2. Synthèse Globale des Écarts
- **Fichiers DAOBS analysés avec observations** : 125
- **Nombre total de pixels-jours comparés** : 416,776
- **Taux d'égalité stricte (`diff == 0.0`)** : **100.0000%** (416,776 / 416,776)
- **Taux d'égalité à précision machine (`diff < 1e-6`)** : **100.0000%** (416,776 / 416,776)
- **Écart absolu maximal** : `0.00000000e+00` $\mathrm{m}^3/\mathrm{m}^3$
- **Écart absolu moyen** : `0.00000000e+00` $\mathrm{m}^3/\mathrm{m}^3$

## 3. Répartition des Écarts
| Seuil d'écart | Nombre de pixels | Pourcentage cumulé |
| :--- | :---: | :---: |
| `diff == 0.0` (stricte égalité) | 416,776 | 100.0000% |
| `diff < 1e-7` | 416,776 | 100.0000% |
| `diff < 1e-6` | 416,776 | 100.0000% |
| `diff >= 1e-6` | 0 | 0.0000% |

## 4. Échantillon Détaillé par Année et Passage
| Fichier DAOBS | Date | Passage | Pixels observés | Égalités exactes | Écart max |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `LISDAOBS_201601010600.a01.d01.1gs4r` | 2016-01-01 | AM | 8837 | 8837 | 0.00e+00 |
| `LISDAOBS_201601010700.a01.d01.1gs4r` | 2016-01-01 | AM | 1973 | 1973 | 0.00e+00 |
| `LISDAOBS_201601011800.a01.d01.1gs4r` | 2016-01-01 | PM | 3390 | 3390 | 0.00e+00 |
| `LISDAOBS_201601011900.a01.d01.1gs4r` | 2016-01-01 | PM | 4 | 4 | 0.00e+00 |
| `LISDAOBS_201601021800.a01.d01.1gs4r` | 2016-01-02 | PM | 1153 | 1153 | 0.00e+00 |
| `LISDAOBS_201601021900.a01.d01.1gs4r` | 2016-01-02 | PM | 1969 | 1969 | 0.00e+00 |
| `LISDAOBS_201601030600.a01.d01.1gs4r` | 2016-01-03 | AM | 11745 | 11745 | 0.00e+00 |
| `LISDAOBS_201601030700.a01.d01.1gs4r` | 2016-01-03 | AM | 120 | 120 | 0.00e+00 |
| `LISDAOBS_201601031800.a01.d01.1gs4r` | 2016-01-03 | PM | 113 | 113 | 0.00e+00 |
| `LISDAOBS_201601040600.a01.d01.1gs4r` | 2016-01-04 | AM | 3151 | 3151 | 0.00e+00 |
| `LISDAOBS_201601040700.a01.d01.1gs4r` | 2016-01-04 | AM | 1970 | 1970 | 0.00e+00 |
| `LISDAOBS_201601050600.a01.d01.1gs4r` | 2016-01-05 | AM | 455 | 455 | 0.00e+00 |
| `LISDAOBS_201601051900.a01.d01.1gs4r` | 2016-01-05 | PM | 58 | 58 | 0.00e+00 |
| `LISDAOBS_201601060600.a01.d01.1gs4r` | 2016-01-06 | AM | 10972 | 10972 | 0.00e+00 |
| `LISDAOBS_201601060700.a01.d01.1gs4r` | 2016-01-06 | AM | 1598 | 1598 | 0.00e+00 |
| `LISDAOBS_201601061800.a01.d01.1gs4r` | 2016-01-06 | PM | 638 | 638 | 0.00e+00 |
| `LISDAOBS_201601070700.a01.d01.1gs4r` | 2016-01-07 | AM | 354 | 354 | 0.00e+00 |
| `LISDAOBS_201601071800.a01.d01.1gs4r` | 2016-01-07 | PM | 6484 | 6484 | 0.00e+00 |
| `LISDAOBS_201601071900.a01.d01.1gs4r` | 2016-01-07 | PM | 1448 | 1448 | 0.00e+00 |
| `LISDAOBS_201601080600.a01.d01.1gs4r` | 2016-01-08 | AM | 6067 | 6067 | 0.00e+00 |
| `LISDAOBS_201601090600.a01.d01.1gs4r` | 2016-01-09 | AM | 7230 | 7230 | 0.00e+00 |
| `LISDAOBS_201601090700.a01.d01.1gs4r` | 2016-01-09 | AM | 1580 | 1580 | 0.00e+00 |
| `LISDAOBS_201601091800.a01.d01.1gs4r` | 2016-01-09 | PM | 3850 | 3850 | 0.00e+00 |
| `LISDAOBS_201601091900.a01.d01.1gs4r` | 2016-01-09 | PM | 34 | 34 | 0.00e+00 |
| `LISDAOBS_201601101800.a01.d01.1gs4r` | 2016-01-10 | PM | 1169 | 1169 | 0.00e+00 |
| `LISDAOBS_201701010700.a01.d01.1gs4r` | 2017-01-01 | AM | 80 | 80 | 0.00e+00 |
| `LISDAOBS_201701011800.a01.d01.1gs4r` | 2017-01-01 | PM | 6107 | 6107 | 0.00e+00 |
| `LISDAOBS_201701011900.a01.d01.1gs4r` | 2017-01-01 | PM | 1785 | 1785 | 0.00e+00 |
| `LISDAOBS_201701020600.a01.d01.1gs4r` | 2017-01-02 | AM | 6589 | 6589 | 0.00e+00 |
| `LISDAOBS_201701030600.a01.d01.1gs4r` | 2017-01-03 | AM | 7798 | 7798 | 0.00e+00 |
| ... (95 autres fichiers analysés) | ... | ... | ... | ... | ... |

## 5. Conclusion
L'écart entre la chaîne de reprojection Python de SAR et les observations injectées dans LIS est **strictement nul** sur l'intégralité du domaine et des périodes testées (couvrant les passages AM et PM de 2016 à 2020). La reprojection par plus proche voisin reproduit avec une fidélité absolue la logique de `read_NASASMAPsm.F90`.
