import pandas as pd

df = pd.read_csv("data/raw.csv", sep=";", low_memory=False)
print("original shape:", df.shape)

keep = [
    "Molecule ChEMBL ID", "Smiles",
    "Standard Type", "Standard Relation", "Standard Value", "Standard Units",
    "pChEMBL Value", "Data Validity Comment", "Potential Duplicate",
    "Assay ChEMBL ID", "Assay Type", "Target Name", "Target Organism",
    "Document ChEMBL ID", "Document Year",
    "Molecular Weight", "AlogP",
]
df = df[keep].copy()
print("after column selection:", df.shape)

for col in ["Standard Type", "Standard Units", "Standard Relation",
            "Assay Type", "Target Name", "Target Organism",
            "Data Validity Comment", "Potential Duplicate"]:
    print(f"\n--- {col} ---")
    print(df[col].value_counts(dropna=False))

print("\n--- counts ---")
print("missing smiles:", df["Smiles"].isna().sum())
print("missing value:", df["Standard Value"].isna().sum())
print("missing pChEMBL:", df["pChEMBL Value"].isna().sum())
print("duplicate IDs:", df["Molecule ChEMBL ID"].duplicated().sum())
print("duplicate SMILES:", df["Smiles"].duplicated().sum())
print("unique compounds:", df["Molecule ChEMBL ID"].nunique())
print("salt/mixture SMILES:", df["Smiles"].str.contains(r"\.", na=False).sum())