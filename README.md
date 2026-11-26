# Compound Data Curation Pipeline

Turns a messy bioactivity export into a clean, standardized, annotated compound table, an activity-cliff list and a client-ready report. One command, fully logged, reproducible.

**Problem it solves.** A discovery team has a compound file assembled from papers, vendors and assays: salts, duplicate structures, mixed units, censored values (`>`, `<`), conflicting measurements. Before any screening or modelling, that has to become one trustworthy row per compound. This tool does that and records every row it removes and why.

**Case study:** ChEMBL IC50 data for CDK4/cyclin D1 (target CHEMBL1907601): 2,204 raw rows -> 1,909 curated compounds. Data downloaded from ChEMBL on [25th Nov 2026, ChEMBL1907601]. ChEMBL data is distributed under CC BY-SA 3.0; check the current terms before redistributing the raw file.

## Quick start

```bash
conda create -n chem -c conda-forge python=3.12 rdkit pandas numpy pytest
conda activate chem

python pipeline.py data/raw.csv -o out
```

This runs all five stages and writes everything to `out/`, ending with `out/report.md`. Use `--skip-report` to skip the report step.

**Input:** a ChEMBL activity export (CSV, semicolon-separated) for one target, IC50 only. Required columns include `Smiles`, `Standard Value`, `Standard Units`, `Standard Relation`, `Molecule ChEMBL ID`, `Target Name`.

**Reproduce the case study:** in ChEMBL, open target CHEMBL1907601 (CDK4/cyclin D1), go to Activities, filter Standard Type = IC50, leave relation and units unfiltered, download as CSV, save as `data/raw.csv`.

## What it does

| Stage | Steps | Output |
|---|---|---|
| 1. Load and standardize | Hard rejects (missing SMILES, value or units; units other than nM), RDKit cleanup, largest-fragment (salt stripping), neutralization, InChIKey | `stage1_cleaned.csv`, `rejected.csv` |
| 2. Activity handling | IC50 (nM) to pIC50, merge duplicate measurements per structure (median), flag high-variance compounds, keep censored-only compounds with the relation flipped to pIC50 direction | `stage2_compounds.csv` |
| 3. Properties and flags | MW, cLogP, TPSA, HBD, HBA, rotatable bonds, Lipinski violations, PAINS alerts | `stage3_annotated.csv` |
| 4. SAR | Murcko scaffolds, Morgan fingerprints, activity cliffs | `stage4_annotated.csv`, `activity_cliffs.csv` |
| 5. Report | Every number computed from the output files, none typed by hand | `report.md` |

Every stage checks that no rows were lost (kept + rejected equals input).

## Case study results

| | |
|---|---|
| Raw rows | 2,204 |
| Rejected | 69 (33 mass units, 31 missing value, 5 missing SMILES) |
| Rows containing salts or mixtures | 44 |
| ChEMBL IDs -> parent structures | 1,911 -> 1,909 |
| Exact measurements -> compounds | 1,762 -> 1,557 (205 duplicate measurements merged) |
| High-variance compounds (range > 1 pIC50) | 11 |
| Censored-only compounds | 352 (246 `<`, 105 `>`, 1 `<=` in pIC50 direction) |
| Compounds with 0 Lipinski violations | 1,341 of 1,909 |
| PAINS alerts | 168 (108 are `anil_di_alk_A`, a common kinase-inhibitor motif) |
| Murcko scaffolds | 848 (599 singletons) |
| Activity cliffs (Tanimoto >= 0.8, delta pIC50 >= 2.0) | 20 pairs, 28 compounds; 9 within one scaffold, 2 stereoisomer pairs |

Largest cliffs:

| A | B | pIC50 A | pIC50 B | Tanimoto |
|---|---|---:|---:|---:|
| CHEMBL280867 | CHEMBL185525 | 4.46 | 8.19 | 0.86 |
| CHEMBL418203 | CHEMBL9130 | 7.05 | 3.60 | 0.86 |

Full report: [`out/report.md`](out/report.md).

## Method notes

- **pIC50** = 9 - log10(IC50 in nM). Validated against ChEMBL's own `pChEMBL Value` on 1,728 exact rows: maximum difference 0.0055, which is ChEMBL's rounding to two decimals.
- **Median, not mean,** when merging duplicates, so a single bad measurement can't drag the value. A compound is flagged `high_variance` when its measurements span more than 1.0 pIC50 (10-fold).
- **Censored values** (`>`, `<`) are not treated as exact. IC50 `> X` means pIC50 `< X`, so the relation is flipped. Compounds that also have exact measurements use only those.
- **Mass-unit rows** (ug/mL) are rejected, not converted: conversion needs molecular weight and is unreliable with salt forms.
- **Cliff definition:** Morgan fingerprint (radius 2, 2048 bits), Tanimoto >= 0.8 and delta pIC50 >= 2.0, on exact, non-high-variance compounds only. Two flags separate cliffs inside one scaffold (`same_scaffold`) and stereoisomer pairs (`same_connectivity`, first block of the InChIKey).
- All thresholds live in `curation/config.py`. The report reads the same constants, so it can't state a different value from the one used.

## Limitations

- Assumes IC50 data for a **single target**. It does not filter by `Standard Type`; the report warns if more than one target is present.
- Standardization is rule-based (RDKit); it does not canonicalize tautomers.
- 1,434 compounds have a single measurement, so the variance check cannot flag them.
- Potencies are medians across assays and papers. Some cliffs may reflect assay differences, not chemistry; assay and document IDs are not carried into the compound table.
- Morgan fingerprints ignore stereochemistry, so stereoisomers score 1.0.
- PAINS alerts and Lipinski violations are flags for review, not filters. Kinase inhibitors often break rule-of-five limits legitimately.
- Only the ChEMBL export format is supported as input.

## Tests

```bash
python -m pytest -q
```

Ten unit tests on small synthetic inputs, with no dependency on the case-study data: salt stripping, charge neutralization, invalid SMILES, mass-unit rejection, row conservation, pIC50 conversion, duplicate merging, relation flipping, Lipinski counting and cliff detection.

## Project layout

```
pipeline.py            command-line entry point
curation/
  config.py            shared thresholds
  load.py              load_raw, hard_rejects
  standardize.py       standardize_smiles
  stage1.py ... stage4.py
  report.py            collect_stats, render_report, write_report
tests/test_curation.py
data/raw.csv           case-study input (ChEMBL export)
out/                   case-study outputs
```

## Tools

Python 3.12, RDKit, pandas, NumPy, pytest.

## License

Code: MIT. Data in `data/` and derived outputs in `out/` come from ChEMBL and follow its licence (CC BY-SA 3.0); please attribute ChEMBL.
