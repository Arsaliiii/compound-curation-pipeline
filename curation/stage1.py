import os
import pandas as pd

from curation.load import hard_rejects, load_raw
from curation.standardize import standardize_smiles


def run_stage1(input_path: str, out_dir: str = "out") -> tuple[pd.DataFrame, pd.DataFrame]:
    # 1. Load raw data
    df = load_raw(input_path)
    n_input = len(df)

    # 2. Apply hard rejects
    kept, rejected = hard_rejects(df)

    # 3. Standardize SMILES; failures go to rejected as "invalid_smiles"
    smiles_results = kept["Smiles"].apply(standardize_smiles)
    invalid_mask = smiles_results.isna()

    if invalid_mask.any():
        invalid_rows = kept[invalid_mask].copy()
        invalid_rows["reject_reason"] = "invalid_smiles"
        rejected = pd.concat([rejected, invalid_rows], ignore_index=True)

    cleaned = kept[~invalid_mask].copy()

    # 4. Add new columns
    valid_tuples = smiles_results[~invalid_mask]
    cleaned["parent_smiles"] = [item[0] for item in valid_tuples]
    cleaned["parent_inchikey"] = [item[1] for item in valid_tuples]
    cleaned["had_salt"] = [item[2] for item in valid_tuples]
    cleaned["censored"] = cleaned["Standard Relation"] != "="

    # 5. Write outputs
    os.makedirs(out_dir, exist_ok=True)
    cleaned.to_csv(os.path.join(out_dir, "stage1_cleaned.csv"), index=False)
    rejected.to_csv(os.path.join(out_dir, "rejected.csv"), index=False)

    # 6. Summary and row-count check
    total = len(cleaned) + len(rejected)
    print(f"Input rows:    {n_input}")
    print(f"Cleaned rows:  {len(cleaned)}")
    print(f"Rejected rows: {len(rejected)}")
    print("Reject breakdown:")
    print(rejected["reject_reason"].value_counts(dropna=False).to_string())
    print(f"Had salt:      {cleaned['had_salt'].sum()}")
    print(f"Censored:      {cleaned['censored'].sum()}")

    assert total == n_input, f"Row count mismatch: {n_input} in, {total} out"

    return cleaned, rejected


if __name__ == "__main__":
    run_stage1("data/raw.csv", "out")