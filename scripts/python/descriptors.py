import pubchempy as pcp
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
    print(f"Ring Count: {rings}")
