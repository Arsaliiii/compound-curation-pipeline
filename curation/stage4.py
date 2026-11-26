import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from curation.config import (
    DELTA_THRESHOLD,
    FINGERPRINT_BITS,
    FINGERPRINT_RADIUS,
    SIM_THRESHOLD,
)

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


fp_gen = rdFingerprintGenerator.GetMorganGenerator(
    radius=FINGERPRINT_RADIUS,
    fpSize=FINGERPRINT_BITS,
)

def _murcko_smiles(smiles: str) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def add_scaffolds(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["scaffold"] = out["parent_smiles"].apply(_murcko_smiles)
    return out

def similar_pairs(df: pd.DataFrame, threshold: float = 0.7) -> pd.DataFrame:
    usable = df[(df["activity_type"] == "exact") & (~df["high_variance"])]

    fps = []
    valid_indices = []
    for idx, smiles in usable["parent_smiles"].items():
        if not isinstance(smiles, str) or not smiles.strip():
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            continue
        fps.append(_MORGAN_GEN.GetFingerprint(mol))
        valid_indices.append(idx)


    pairs: list[dict[str, int | float]] = []
    n = len(fps)
    for i in range(n - 1):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :])
        for offset, sim in enumerate(sims):
            if sim >= threshold:
                pairs.append(
                    {
                        "idx_a": valid_indices[i],
                        "idx_b": valid_indices[i + 1 + offset],
                        "similarity": float(sim),
                    }
                )

    return pd.DataFrame(pairs, columns=["idx_a", "idx_b", "similarity"])

def find_activity_cliffs(
    df: pd.DataFrame,
    sim_threshold: float = 0.8,
    delta_threshold: float = 2.0,
) -> pd.DataFrame:
    output_cols = [
        "idx_a", "idx_b",
        "chembl_id_a", "chembl_id_b",
        "inchikey_a", "inchikey_b",
        "smiles_a", "smiles_b",
        "pIC50_a", "pIC50_b",
        "similarity", "delta_pIC50",
    ]


    pairs = similar_pairs(df, sim_threshold)
    if pairs.empty:
        return pd.DataFrame(columns=output_cols)

    idx_a = pairs["idx_a"]
    idx_b = pairs["idx_b"]

    p_a = df.loc[idx_a, "pIC50"].to_numpy(dtype=float)
    p_b = df.loc[idx_b, "pIC50"].to_numpy(dtype=float)
    delta = np.abs(p_a - p_b)

    cliff_mask = np.isfinite(delta) & (delta >= delta_threshold)
    ia = idx_a[cliff_mask]
    ib = idx_b[cliff_mask]

    cliffs = pd.DataFrame(
        {
            "idx_a": ia.to_numpy(),
            "idx_b": ib.to_numpy(),
            "chembl_id_a": df.loc[ia, "chembl_id"].to_numpy(),
            "chembl_id_b": df.loc[ib, "chembl_id"].to_numpy(),
            "inchikey_a": df.loc[ia, "parent_inchikey"].to_numpy(),
            "inchikey_b": df.loc[ib, "parent_inchikey"].to_numpy(),
            "smiles_a": df.loc[ia, "parent_smiles"].to_numpy(),
            "smiles_b": df.loc[ib, "parent_smiles"].to_numpy(),
            "pIC50_a": p_a[cliff_mask],
            "pIC50_b": p_b[cliff_mask],
            "similarity": pairs.loc[cliff_mask, "similarity"].to_numpy(dtype=float),
            "delta_pIC50": delta[cliff_mask],
        }
    )

    cliffs = cliffs.sort_values("delta_pIC50", ascending=False).reset_index(drop=True)
    return cliffs[output_cols]

def run_stage4(
    input_path: str = "out/stage3_annotated.csv",
    out_dir: str = "out",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(input_path)
    n_input = len(df)

    stage4_annotated = add_scaffolds(df)

    os.makedirs(out_dir, exist_ok=True)
    stage4_annotated.to_csv(os.path.join(out_dir, "stage4_annotated.csv"), index=False)

    cliffs = find_activity_cliffs(stage4_annotated, SIM_THRESHOLD, DELTA_THRESHOLD)

    cliffs["same_connectivity"] = (
        cliffs["inchikey_a"].astype(str).str[:14] == cliffs["inchikey_b"].astype(str).str[:14]
    )
    scaff_a = stage4_annotated.loc[cliffs["idx_a"], "scaffold"].to_numpy()
    scaff_b = stage4_annotated.loc[cliffs["idx_b"], "scaffold"].to_numpy()
    cliffs["same_scaffold"] = (scaff_a != "") & (scaff_a == scaff_b)

    cliffs.to_csv(os.path.join(out_dir, "activity_cliffs.csv"), index=False)

    scaffold_counts = stage4_annotated["scaffold"].value_counts()
    n_cliff_compounds = pd.concat([cliffs["chembl_id_a"], cliffs["chembl_id_b"]]).nunique()

    print(f"Unique scaffolds: {len(scaffold_counts)} (singletons: {int((scaffold_counts == 1).sum())})")
    print(f"Cliffs found: {len(cliffs)} (unique compounds involved: {n_cliff_compounds})")
    print(f"same_connectivity == True: {int(cliffs['same_connectivity'].sum())}")
    print(f"same_scaffold == True: {int(cliffs['same_scaffold'].sum())}")

    assert len(stage4_annotated) == n_input, (
        f"Expected {n_input} rows in stage4_annotated, got {len(stage4_annotated)}"
    )

    return stage4_annotated, cliffs

if __name__ == "__main__":
    run_stage4("out/stage3_annotated.csv", "out")