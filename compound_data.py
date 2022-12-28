from rdkit import Chem
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
import requests


class Compound:
    def __init__(self):
        self.compound_img = None
        self.name = ""
        self.smiles = ""
        self.inchikey = ""
        self.num_of_atoms = ""
        self.mol_weight = ""
        self.num_of_bonds = ""
        self.num_of_rotatable_bonds = ""

    def smiles_to_iupac(self, smiles):
        name = ""
        try:
            CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
            rep = "iupac_name"
            url = CACTUS.format(smiles, rep)
            response = requests.get(url)
            response.raise_for_status()
            name = response.text
        except:
            name = ""
        return name

    def create_compound_from_smile(self, input):
        chem_compound = Chem.MolFromSmiles(input)
        self.input_data_from_compound(chem_compound)

    def parse_compound_from_link(self, url):
        chem_mol = requests.get(url).text
        chem_compound = Chem.MolFromMolBlock(chem_mol)
        self.input_data_from_compound(chem_compound)

    def input_data_from_compound(self, chem_compound):
        self.compound_img = Draw.MolToImage(chem_compound)
        self.smiles = Chem.MolToSmiles(chem_compound)
        self.name = self.smiles_to_iupac(self.smiles)
        self.inchikey = Chem.MolToInchiKey(chem_compound)
        self.num_of_atoms = chem_compound.GetNumAtoms()
        # print([atom.GetSymbol() for atom in chem_compound.GetAtoms()])
        self.mol_weight = Descriptors.MolWt(chem_compound)
        self.num_of_bonds = chem_compound.GetNumBonds()
        self.num_of_rotatable_bonds = Descriptors.NumRotatableBonds(chem_compound)
