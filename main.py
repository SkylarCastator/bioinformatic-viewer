import streamlit as st
#import sequence_analyzer as covid

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import requests

class ui_instance:
    def __init__(self):
        self.input_chemical_path = "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol"

        st.write("""
        # Bio-Informatic Reader
        Testing a user interface to view information of medicine and diseases.
        """)
        self.tools_sidebar()
        self.analyze_chemical_compond()

    def analyze_chemical_compond(self):
        chem_mol = requests.get(self.input_chemical_path).text
        chem_compound = Chem.MolFromMolBlock(chem_mol)
        im = Draw.MolToImage(chem_compound)
        print(chem_compound.GetNumAtoms())
        # print([atom.GetSymbol() for atom in chem_compound.GetAtom()])
        #print(Descriptors.MolWt(chem_compound))
        #print(Descriptors.NumRotatableBonds(chem_compound))
        st.image(im, caption='Analysis of Chemical')

    def tools_sidebar(self):
        st.sidebar.write("""
        # Toolbar
        Adjust these fields to learn more!
        """)
        chemical_path = st.sidebar.text_input(
            "Chemical Compound Path",
            "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol")
        self.input_chemical_path = chemical_path
        if st.sidebar.button('Search'):
            self.analyze_chemical_compond()




#Id NC_045512.2
#https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=graph
#https://www.rcsb.org/structure/6ZCO
#nv.demo()
#https://william-dawson.github.io/using-py3dmol.html
#https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753
#sequence_file_path = "../data_samples/Covid_sequence-NC_045512.fasta"
#sequence_data = covid.SequenceAnalyzer(sequence_file_path)
ui_instance = ui_instance()