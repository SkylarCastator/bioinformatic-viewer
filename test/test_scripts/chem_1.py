import streamlit as st # https://streamlit.io/
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

chemical_smile = 'CN1CCC[C@H]1c2cccnc2'
Chem.MolFromSmiles(chemical_smile)
#https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
compound_smiles = 'c1cc(C(=O)O)c(OC(=O)C)cc1'
m = Chem.MolFromSmiles(compound_smiles)
im=Draw.MolToImage(m)
im.show()

#st.image(im)

import requests
chem_link = "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol"
chem_mol = requests.get(chem_link).text
chem_compound = Chem.MolFromMolBlock(chem_mol)
im=Draw.MolToImage(chem_compound)
im.show()
print(chem_compound.GetNumAtoms())
print([atom.GetSymbol() for atom in chem_compound.GetAtom()])

from rdkit.Chem import Descriptors
print(Descriptors.MolWt(chem_compound))
print(Descriptors.NumRotatableBonds(chem_compound))

