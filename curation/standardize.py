from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

def standardize_smiles(smiles: str) -> tuple[str, str, bool] | None:
    # returns (parent_smiles, inchikey, had_salt), or None if the SMILES can't be parsed
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    had_salt = "." in smiles or len(Chem.GetMolFrags(mol)) > 1

    mol = rdMolStandardize.Cleanup(mol)                              # added
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    mol = rdMolStandardize.Uncharger().uncharge(mol)                 # added

    parent_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchikey = MolToInchiKey(mol)

    return (parent_smiles, inchikey, had_salt)