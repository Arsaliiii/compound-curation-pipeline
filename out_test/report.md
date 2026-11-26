# Bioactivity Curation & SAR Report

## 1. Summary

- Target: Cyclin-dependent kinase 4/cyclin D1 (IC50, binding assays).
- 2,204 input rows: 2,135 cleaned, 69 rejected.
- 1,909 compounds: 1,557 with exact potency, 352 censored-only.

## 2. Rejected rows

| Reason | Count |
|---|---|
| non_molar_units | 33 |
| missing_value | 31 |
| missing_smiles | 5 |

## 3. Standardization

- 44 rows contained salts or mixtures (largest fragment kept).
- 1,911 ChEMBL IDs map to 1,909 parent structures.

## 4. Activity data

- 1,762 exact measurements merged into 1,557 compounds (205 duplicate measurements merged).
- 11 compounds flagged: measurements span more than 1.0 pIC50.

| Censored-only relation (pIC50) | Count |
|---|---|
| < | 246 |
| > | 105 |
| <= | 1 |

## 5. Properties and flags

- Mean MW 440.5, mean cLogP 3.69.
- 168 compounds match a PAINS pattern.

| Lipinski violations | Compounds |
|---|---|
| 0 | 1,341 |
| 1 | 365 |
| 2 | 183 |
| 3 | 19 |
| 4 | 1 |

| Top PAINS patterns | Compounds |
|---|---|
| anil_di_alk_A(478) | 108 |
| indol_3yl_alk(461) | 24 |
| catechol_A(92) | 8 |

## 6. SAR and activity cliffs

- 848 Murcko scaffolds (599 singletons).
- 20 cliff pairs involving 28 compounds (9 within one scaffold, 2 stereoisomer pairs).

| Top scaffolds | Compounds |
|---|---|
| `O=C1NC(=O)c2c1c1c(ccc3cc[nH]c31)c1[nH]c3ccccc3c21` | 48 |
| `O=c1[nH]cnc2c1cnn2-c1ccccc1` | 47 |
| `O=c1[nH]c(Cc2ccccc2)nc2c1cnn2-c1ccccc1` | 28 |
| `O=C1c2ccccc2-c2n[nH]c(-c3ccccc3)c21` | 26 |
| `O=C1NCC2(CCCCC2)n2c1cc1cnc(Nc3ccc(N4CCNCC4)cn3)nc12` | 25 |

| A | B | pIC50 A | pIC50 B | Similarity | Delta | Same scaffold |
|---|---|---|---|---|---|---|
| CHEMBL280867 | CHEMBL185525 | 4.46 | 8.19 | 0.86 | 3.73 | True |
| CHEMBL418203 | CHEMBL9130 | 7.05 | 3.60 | 0.86 | 3.44 | False |
| CHEMBL300490 | CHEMBL9130 | 7.00 | 3.60 | 0.86 | 3.40 | False |
| CHEMBL180431 | CHEMBL361769 | 9.00 | 5.76 | 0.85 | 3.24 | False |
| CHEMBL418203 | CHEMBL8993 | 7.05 | 3.92 | 0.86 | 3.12 | False |

## 7. Parameters

- High-variance range: more than 1.0 pIC50
- Fingerprint: Morgan, radius 2, 2048 bits
- Activity cliff: Tanimoto at least 0.8 and delta pIC50 at least 2.0

## 8. Limitations

- 1,434 compounds have one measurement, so the variance check cannot flag them.
- 352 censored-only compounds are excluded from variance and cliff analysis.
- 11 high-variance compounds are also excluded from cliff analysis.
- Potencies are medians across assays, so some cliffs may reflect assay differences, not chemistry.
- Morgan fingerprints ignore stereochemistry, so stereoisomers can look identical.
- anil_di_alk_A(478) accounts for 108 of 168 PAINS flags; flags are alerts, not proof.
- 33 rows reported in mass units (ug/mL) were rejected, not converted.
