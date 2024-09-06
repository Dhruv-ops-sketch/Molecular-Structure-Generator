from rdkit import Chem
from rdkit.Chem import Draw
from pubchempy import get_compounds
import os

def compound_to_smiles(compound_name):
    compounds = get_compounds(compound_name, 'name')
    if compounds:
        smiles = compounds[0].canonical_smiles
        return smiles
    else:
        return "Compound not found"

def draw_molecule(smiles, compound_name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Failed to generate molecular structure")
    else:
        folder_path = os.path.join(os.path.dirname(__file__), 'Generated Structures')
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        image_path = os.path.join(folder_path, f'{compound_name}.png')
        Draw.MolToFile(mol, image_path)

compound_name = input("Enter the compound name: ")
smiles = compound_to_smiles(compound_name)
draw_molecule(smiles, compound_name)
