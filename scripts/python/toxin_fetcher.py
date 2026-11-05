import pubchempy as pcp

print("Fetching data from PubChem...")
# get_compounds returns a list, so we take the first result [0]
results = pcp.get_compounds('tetrodotoxin', 'name')
toxin = results[0]

print(f"Name: Tetrodotoxin")
print(f"Formula: {toxin.molecular_formula}")
print(f"Molecular Weight: {toxin.molecular_weight}")
print(f"SMILES: {toxin.isomeric_smiles}")
