from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from pubchempy import get_compounds
import py3Dmol

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
        img = Draw.MolToImage(mol)
        return img

def generate_3d_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Failed to generate 3D structure from SMILES")
        return None
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
    AllChem.UFFOptimizeMolecule(mol)  # Optimize the structure
    return mol

def visualize_3d_structure(mol):
    block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(block, 'mol')  # Add molecule to viewer
    viewer.setStyle({'stick': {}})  # Set the style to "stick" representation
    viewer.zoomTo()  # Zoom to fit the molecule
    return viewer.show()

# User interaction in Jupyter
compound_name = input("Enter the compound name: ")
smiles = compound_to_smiles(compound_name)
if smiles != "Compound not found":
    img = draw_molecule(smiles, compound_name)
    img.show()  # Show 2D structure in Jupyter

    # Generate and visualize 3D structure
    mol_3d = generate_3d_structure(smiles)
    if mol_3d:
        visualize_3d_structure(mol_3d)
else:
    print(smiles)
