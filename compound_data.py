from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import requests


class Compound:
    def __init__(self):
        self.compound_img = None
        self.num_of_atoms = ""
        self.mol_weight = ""
        self.num_of_rotatable_bonds = ""

    def parse_compound_from_link(self, url):
        chem_mol = requests.get(url).text
        chem_compound = Chem.MolFromMolBlock(chem_mol)
        self.compound_img = Draw.MolToImage(chem_compound)
        self.num_of_atoms = chem_compound.GetNumAtoms()
        # print([atom.GetSymbol() for atom in chem_compound.GetAtom()])
        self.mol_weight = Descriptors.MolWt(chem_compound)
        self.num_of_rotatable_bonds = Descriptors.NumRotatableBonds(chem_compound)