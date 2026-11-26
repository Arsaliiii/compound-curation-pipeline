import pandas as pd
import pytest

from curation.load import hard_rejects
from curation.stage2 import FLIP, add_pic50, merge_exact
from curation.stage3 import compute_properties, lipinski_violations
from curation.stage4 import find_activity_cliffs
from curation.standardize import standardize_smiles


# 1. Salt stripped
def test_salt_is_stripped():
    parent, _, had_salt = standardize_smiles("CC(=O)O.[Na]")
    assert parent == "CC(=O)O"
    assert had_salt is True


# 2. Charge neutralized & same parent InChIKey as test 1
def test_charge_neutralized_same_inchikey():
    _, key_salt, _ = standardize_smiles("CC(=O)O.[Na]")
    parent_zwit_or_ion, key_ion, _ = standardize_smiles("CC(=O)[O-]")
    assert parent_zwit_or_ion == "CC(=O)O"
    assert key_ion == key_salt


# 3. Invalid SMILES returns None
def test_invalid_smiles_returns_none():
    assert standardize_smiles("not_a_smiles") is None


# 4. Mass units rejected with reason non_molar_units
def test_mass_units_rejected():
    df = pd.DataFrame(
        {
            "Smiles": ["CC(=O)Oc1ccccc1C(=O)O"],
            "Standard Value": [10.0],
            "Standard Units": ["ug.mL-1"],
            "Standard Relation": ["="],
            "Molecule ChEMBL ID": ["CHEMBL_MASS_1"],
        }
    )
    kept, rejected = hard_rejects(df)
    assert len(kept) == 0
    assert len(rejected) == 1
    assert rejected.iloc[0]["reject_reason"] == "non_molar_units"


# 5. Row conservation through hard_rejects
def test_row_conservation_hard_rejects():
    df = pd.DataFrame(
        {
            "Smiles": [
                "CC(=O)Oc1ccccc1C(=O)O",
                "CCO",
                "c1ccccc1",
                "CCN(CC)CC",
            ],
            "Standard Value": [100.0, 5.0, 20.0, 50.0],
            "Standard Units": ["nM", "ug.mL-1", "nM", "uM"],
            "Standard Relation": ["=", "=", "=", "="],
            "Molecule ChEMBL ID": ["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4"],
        }
    )
    kept, rejected = hard_rejects(df)
    assert len(kept) + len(rejected) == len(df)


# 6. pIC50 conversion: 100 nM -> 7.0, 1 nM -> 9.0
def test_pic50_conversion():
    df = pd.DataFrame(
        {
            "Standard Value": [100.0, 1.0],
            "Standard Units": ["nM", "nM"],
        }
    )
    out = add_pic50(df)
    assert out["pIC50"].iloc[0] == pytest.approx(7.0)
    assert out["pIC50"].iloc[1] == pytest.approx(9.0)


# 7. Merge exact duplicates: median, n_measurements, and high_variance flag
def test_merge_exact_high_variance():
    _, key, _ = standardize_smiles("CC(=O)Oc1ccccc1C(=O)O")
    df = pd.DataFrame(
        {
            "censored": [False, False],
            "pIC50": [6.0, 8.0],
            "parent_inchikey": [key, key],
            "parent_smiles": ["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
            "Molecule ChEMBL ID": ["CHEMBL_A", "CHEMBL_A"],
        }
    )
    merged = merge_exact(df)
    assert len(merged) == 1
    row = merged.iloc[0]
    assert int(row["n_measurements"]) == 2
    assert row["pIC50_median"] == pytest.approx(7.0)
    assert bool(row["high_variance"]) is True


# 8. Relation flip for censored-only rows (> and < become < and >)
def test_relation_flip():
    relations = pd.Series([">", "<"])
    flipped = relations.map(FLIP).tolist()
    assert flipped == ["<", ">"]


# 9. Lipinski violations: MW 600 and HBD 6 -> 2 violations
def test_lipinski_violations_count():
    props = {
        "MW": 600.0,
        "cLogP": 2.5,
        "TPSA": 80.0,
        "HBD": 6,
        "HBA": 4,
        "RotB": 5,
    }
    assert lipinski_violations(props) == 2

    # Sanity check that compute_properties returns expected keys on a valid SMILES
    aspirin_props = compute_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin_props is not None
    assert lipinski_violations(aspirin_props) == 0


# 10. Activity cliff detection on high-similarity pair with delta >= 2.0
def test_find_activity_cliffs_single_cliff():
    s_l, k_l, _ = standardize_smiles("C[C@H](N)C(=O)O")
    s_d, k_d, _ = standardize_smiles("C[C@@H](N)C(=O)O")
    s_far, k_far, _ = standardize_smiles("CCCCCCCCCCCCCCCC")

    df = pd.DataFrame(
        {
            "activity_type": ["exact", "exact", "exact"],
            "high_variance": [False, False, False],
            "parent_smiles": [s_l, s_d, s_far],
            "parent_inchikey": [k_l, k_d, k_far],
            "pIC50": [5.0, 8.0, 6.0],
            "chembl_id": ["CHEMBL_1", "CHEMBL_2", "CHEMBL_3"],
        }
    )

    cliffs = find_activity_cliffs(df)
    assert len(cliffs) == 1
    row = cliffs.iloc[0]
    assert {row["chembl_id_a"], row["chembl_id_b"]} == {"CHEMBL_1", "CHEMBL_2"}
    assert row["delta_pIC50"] == pytest.approx(3.0)
    assert row["similarity"] == pytest.approx(1.0)