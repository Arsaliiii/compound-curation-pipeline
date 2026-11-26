import os

import numpy as np
import pandas as pd

from curation.config import HIGH_VARIANCE_RANGE

def add_pic50(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    values = pd.to_numeric(out["Standard Value"], errors="coerce")
    valid_mask = values.notna() & (values > 0)

    out["pIC50"] = np.nan
    out.loc[valid_mask, "pIC50"] = 9.0 - np.log10(values[valid_mask].to_numpy(dtype=float))

    return out

def merge_exact(df: pd.DataFrame) -> pd.DataFrame:
    exact = df[(~df["censored"]) & df["pIC50"].notna()].copy()

    merged = (
        exact.groupby("parent_inchikey", as_index=False)
        .agg(
            n_measurements=("pIC50", "count"),
            pIC50_median=("pIC50", "median"),
            pIC50_min=("pIC50", "min"),
            pIC50_max=("pIC50", "max"),
            parent_smiles=("parent_smiles", "first"),
            chembl_id=("Molecule ChEMBL ID", "first"),
        )
    )

    merged["pIC50_range"] = merged["pIC50_max"] - merged["pIC50_min"]
    merged["high_variance"] = merged["pIC50_range"] > HIGH_VARIANCE_RANGE

    return merged

FLIP = {">": "<", ">=": "<=", "<": ">"}


def merge_censored(df: pd.DataFrame, merged: pd.DataFrame) -> pd.DataFrame:
    censored_only = df[
        df["censored"] & ~df["parent_inchikey"].isin(merged["parent_inchikey"])
    ].copy()

    censored_only["pIC50_relation"] = censored_only["Standard Relation"].map(FLIP)

    censored_merged = (
        censored_only.groupby("parent_inchikey", as_index=False)
        .agg(
            parent_smiles=("parent_smiles", "first"),
            chembl_id=("Molecule ChEMBL ID", "first"),
            pIC50_bound=("pIC50", "first"),
            pIC50_relation=("pIC50_relation", "first"),
            n_measurements=("pIC50", "count"),
            n_relations=("pIC50_relation", "nunique"),
        )
    )

    columns = [
        "parent_inchikey",
        "parent_smiles",
        "chembl_id",
        "pIC50_bound",
        "pIC50_relation",
        "n_measurements",
        "n_relations",
    ]
    return censored_merged[columns]

def run_stage2(input_path: str = "out/stage1_cleaned.csv", out_dir: str = "out") -> pd.DataFrame:
    df = pd.read_csv(input_path)
    df = add_pic50(df)

    merged = merge_exact(df)
    censored = merge_censored(df, merged)

    exact_part = merged.assign(
        activity_type="exact",
        pIC50=merged["pIC50_median"],
        pIC50_relation="=",
    )

    censored_part = censored.assign(
        activity_type="censored_only",
        pIC50=censored["pIC50_bound"],
        high_variance=False,
    )

    result = pd.concat([exact_part, censored_part], ignore_index=True)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "stage2_compounds.csv")
    result.to_csv(out_path, index=False)

    print(f"Total unique compounds: {len(result)}")
    print(result["activity_type"].value_counts(dropna=False).to_string())
    print(f"High-variance exact compounds: {int(result['high_variance'].sum())}")

    assert len(result) == df["parent_inchikey"].nunique(), (
        f"Expected {df['parent_inchikey'].nunique()} compounds, got {len(result)}"
    )

    return result

if __name__ == "__main__":
    run_stage2("out/stage1_cleaned.csv", "out")
