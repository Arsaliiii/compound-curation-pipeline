import pandas as pd

KEEP_COLUMNS = [
    "Molecule ChEMBL ID", "Smiles",
    "Standard Type", "Standard Relation", "Standard Value", "Standard Units",
    "pChEMBL Value", "Data Validity Comment", "Potential Duplicate",
    "Assay ChEMBL ID", "Assay Type", "Target Name", "Target Organism",
    "Document ChEMBL ID", "Document Year",
    "Molecular Weight", "AlogP",
]

def load_raw(path: str, keep: list[str] | None = None) -> pd.DataFrame:
    cols = keep if keep is not None else KEEP_COLUMNS
    df = pd.read_csv(path, sep=";", usecols=lambda c: c in cols)
    df = df[[c for c in cols if c in df.columns]].copy()

    if "Standard Relation" in df.columns:
        df["Standard Relation"] = (
            df["Standard Relation"].astype(str).str.strip("'\"")
        )

    return df


def hard_rejects(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = df.copy()
    out["reject_reason"] = None

    empty_reason = out["reject_reason"].isna()

    missing_smiles = empty_reason & out["Smiles"].isna()
    out.loc[missing_smiles, "reject_reason"] = "missing_smiles"

    empty_reason = out["reject_reason"].isna()
    missing_value = empty_reason & out["Standard Value"].isna()
    out.loc[missing_value, "reject_reason"] = "missing_value"

    empty_reason = out["reject_reason"].isna()
    missing_units = empty_reason & out["Standard Units"].isna()
    out.loc[missing_units, "reject_reason"] = "missing_units"

    empty_reason = out["reject_reason"].isna()
    non_molar = empty_reason & (out["Standard Units"] != "nM")
    out.loc[non_molar, "reject_reason"] = "non_molar_units"

    rejected = out[out["reject_reason"].notna()].copy()
    kept = out[out["reject_reason"].isna()].drop(columns=["reject_reason"]).copy()

    return kept, rejected