import os
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

# Ensure the images directory exists so the script doesn't crash
os.makedirs('images', exist_ok=True)

toxin_names = ['tetrodotoxin', 'batrachotoxin']

for name in toxin_names:
    print(f"Drawing {name}...")
    
    # Fetch SMILES from PubChem
    results = pcp.get_compounds(name, 'name')
    smiles = results[0].isomeric_smiles
    
    # Convert text SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    # Save as a high-resolution PNG file
    filename = f"images/{name}_2d.png"
    Draw.MolToFile(mol, filename, size=(600, 600))
    
    print(f"Saved {filename}")
