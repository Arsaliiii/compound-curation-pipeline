import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

_params = FilterCatalogParams()
_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
PAINS_CATALOG = FilterCatalog(_params)


def compute_properties(smiles: str) -> dict[str, float | int] | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
    }


def lipinski_violations(props: dict[str, float | int] | None) -> int:
    if not props:
        return 0

    violations = 0
    if props["MW"] > 500:
        violations += 1
    if props["cLogP"] > 5:
        violations += 1
    if props["HBD"] > 5:
        violations += 1
    if props["HBA"] > 10:
        violations += 1

    return violations


def pains_match(smiles: str) -> str | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    entry = PAINS_CATALOG.GetFirstMatch(mol)
    if entry is None:
        return None

    return entry.GetDescription()

def run_stage3(
    input_path: str = "out/stage2_compounds.csv",
    out_dir: str = "out",
) -> pd.DataFrame:
    df = pd.read_csv(input_path)

    props_list = df["parent_smiles"].apply(compute_properties).tolist()
    props_df = pd.DataFrame(props_list, index=df.index)

    out = pd.concat([df, props_df], axis=1)
    out["lipinski_violations"] = [lipinski_violations(p) for p in props_list]
    out["pains_name"] = out["parent_smiles"].apply(pains_match)
    out["pains_flag"] = out["pains_name"].notna() & (out["pains_name"].astype(str).str.len() > 0)

    os.makedirs(out_dir, exist_ok=True)
    out.to_csv(os.path.join(out_dir, "stage3_annotated.csv"), index=False)

    print(f"Annotated compounds: {len(out)}")
    print("Lipinski violations breakdown:")
    print(out["lipinski_violations"].value_counts().sort_index().to_string())
    print(f"PAINS flagged compounds: {int(out['pains_flag'].sum())}")
    print("Top PAINS patterns:")
    print(out["pains_name"].value_counts().head(3).to_string())
    print(out[["MW", "cLogP"]].describe().loc[["mean", "min", "max"]].to_string())
    print("Missing property values:", int(out[["MW","cLogP","TPSA","HBD","HBA","RotB"]].isna().sum().sum()))

    n_input = len(df)
    ...
    assert len(out) == n_input, f"Expected {n_input} rows, got {len(out)}"

    return out


if __name__ == "__main__":
    run_stage3("out/stage2_compounds.csv", "out")