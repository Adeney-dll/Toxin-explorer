"""import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# The two toxins to compare
toxin_names = ['tetrodotoxin', 'batrachotoxin']

for name in toxin_names:
    print(f"\n--- {name.capitalize()} ---")
    
    # Fetch SMILES from PubChem
    results = pcp.get_compounds(name, 'name')
    smiles = results[0].isomeric_smiles
    
    # Convert text SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    # Calculate descriptors
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Print results
    print(f"Molecular Weight: {mw:.2f}")
    print(f"LogP (Lipophilicity): {logp:.2f}")
    print(f"H-Bond Donors: {hbd}")
    print(f"H-Bond Acceptors: {hba}")
    print(f"Polar Surface Area: {tpsa:.2f}")
    print(f"Ring Count: {rings}")"""


from pathlib import Path

import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parent.parent.parent
INPUT_FILE = ROOT / "data" / "toxins.txt"
OUTPUT_FILE = ROOT / "data" / "toxin_descriptors.csv"


def read_names(path):
    """One toxin name per line; blank lines are ignored."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return [line.strip() for line in lines if line.strip()]


def fetch_compound(name):
    """Return (cid, formula, smiles) for the first PubChem match, or None."""
    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        return None
    c = compounds[0]
    # Newer PubChemPy versions use .smiles; older ones use .isomeric_smiles
    smiles = getattr(c, "smiles", None) or c.isomeric_smiles
    return c.cid, c.molecular_formula, smiles


def calculate_descriptors(mol):
    return {
        "mol_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "h_donors": rdMolDescriptors.CalcNumHBD(mol),
        "h_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
    }


def main():
    rows = []
    failed = []

    for name in read_names(INPUT_FILE):
        print(f"Processing {name}...")
        try:
            record = fetch_compound(name)
            if record is None:
                raise ValueError("no PubChem match")
            cid, formula, smiles = record

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("RDKit could not parse the SMILES")

            row = {"name": name, "cid": cid, "formula": formula}
            row.update(calculate_descriptors(mol))
            rows.append(row)
        except Exception as err:
            print(f"  FAILED: {name} ({err})")
            failed.append(name)

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_FILE, index=False)  # overwritten each run, never appended
    print(f"\nWrote {len(df)} rows to {OUTPUT_FILE}")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
