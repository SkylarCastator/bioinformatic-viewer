import streamlit as st
import pandas as pd
import compound_data
#import sequence_analyzer as covid

class ui_instance:
    def __init__(self):
        self.input_chemical_path = "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol"
        self.chemical_image = None
        self.chemical_title = None
        self.header = st.write("""
        # Bio-Informatic Reader
        Testing a user interface to view information of medicine and diseases.
        """)
        self.tools_sidebar()
        #self.analyze_chemical_compond()

    def analyze_chemical_compond(self, compound_instance):
        left_column, right_column = st.columns(2)
        with left_column:
            self.chemical_image = st.image(compound_instance.compound_img, caption='Analysis of Chemical')
        with right_column:
            self.chemical_title = st.write("## Chemical Name")
            d = {
                'Key': ['Name', 'Smiles', 'Inchikey', 'Num of Atoms', 'Mol Weight', 'Num of Rotatable Bonds'],
                'Value': [
                    compound_instance.name,
                    compound_instance.smiles,
                    compound_instance.inchikey,
                    compound_instance.num_of_atoms,
                    compound_instance.mol_weight,
                    compound_instance.num_of_rotatable_bonds]
            }
            df = pd.DataFrame(data=d)
            st.dataframe(df)

    def tools_sidebar(self):
        st.sidebar.write("""
        # Toolbar
        Adjust these fields to learn more!
        """)
        st.sidebar.write("""
                ## DNA Tools
                """)
        st.sidebar.write("""
                ## Chemical Compound Tools
                """)

        chemical_path = st.sidebar.text_input(
            "Chemical Compound Path",
            "https://go.drugbank.com/structures/small_molecule_drugs/DB14761.mol")
        self.input_chemical_path = chemical_path
        if st.sidebar.button('Search'):
            compound_instance = compound_data.Compound()
            compound_instance.parse_compound_from_link(self.input_chemical_path)
            self.analyze_chemical_compond(compound_instance)

        smile_input = st.sidebar.text_input(
            "Smiles Input",
            'CN1CCC[C@H]1c2cccnc2')
        if st.sidebar.button('Create'):
            compound_instance = compound_data.Compound()
            compound_instance.create_compound_from_smile(smile_input)
            self.analyze_chemical_compond(compound_instance)




#Id NC_045512.2
#https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=graph
#https://www.rcsb.org/structure/6ZCO
#nv.demo()
#https://william-dawson.github.io/using-py3dmol.html
#https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753
#sequence_file_path = "../data_samples/Covid_sequence-NC_045512.fasta"
#sequence_data = covid.SequenceAnalyzer(sequence_file_path)
ui_instance = ui_instance()