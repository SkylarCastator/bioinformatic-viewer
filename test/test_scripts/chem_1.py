from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

#chemical_smile = 'CN1CCC[C@H]1c2cccnc2'
#Chem.MolFromSmiles(chemical_smile)
#https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
#compound_smiles = 'c1cc(C(=O)O)c(OC(=O)C)cc1'
#m = Chem.MolFromSmiles(compound_smiles)
#im=Draw.MolToImage(m)
#im.show()

#st.image(im)

import requests
chem_link = "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol"
chem_mol = requests.get(chem_link).text
print(chem_mol)
chem_compound = Chem.MolFromMolBlock(chem_mol)
im=Draw.MolToImage(chem_compound)
#im.show()
print(chem_compound.GetNumAtoms())
print([atom.GetSymbol() for atom in chem_compound.GetAtoms()])
print(chem_compound.GetNumBonds())
print(chem_compound.GetPropNames())
print(chem_compound.GetNumHeavyAtoms())

print(Descriptors.MolWt(chem_compound))
print(Descriptors.NumRotatableBonds(chem_compound))



CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"

def smiles_to_iupac(smiles):
    rep = "iupac_name"
    url = CACTUS.format(smiles, rep)
    response = requests.get(url)
    response.raise_for_status()
    return response.text


print(smiles_to_iupac('c1ccccc1'))
print(smiles_to_iupac('CC(=O)OC1=CC=CC=C1C(=O)O'))

try:
    smiles = Chem.MolToSmiles(chem_compound)
    print(smiles)
    smiles = smiles.replace('#', '%23')
    #smiles.replace('C#C', 'C%23C')
    print(smiles_to_iupac(smiles))
except:
    pass

#https://xinhaoli74.github.io/posts/2020/04/RDKit-Cheatsheet/
#https://www.rdkit.org/docs/source/rdkit.Chem.PropertyMol.html
#https://github.com/rdkit/rdkit-tutorials/blob/master/notebooks/005_Chemical_space_analysis_and_visualization.ipynb
#https://www.rdkit.org/docs/GettingStartedInPython.html
#https://www.rdkit.org/docs/Cookbook.html