import streamlit as st
import pandas as pd
import compound_data
import dna_data
from collections import Counter

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

    def analyze_dna_sequence(self, dna_sequence):
        st.subheader('Nitrogenous Bases')
        n_count = Counter(dna_sequence.dna)
        df = pd.DataFrame.from_records(list(dict(n_count).items()), columns=['bases', 'count'])
        st.bar_chart(df)

        st.subheader('Proteins')
        n_count = Counter(dna_sequence.proteins)
        df = pd.DataFrame.from_records(list(dict(n_count).items()), columns=['proteins', 'count'])
        st.bar_chart(df)


    def tools_sidebar(self):
        st.sidebar.write("""
        # Toolbar
        Adjust these fields to learn more!
        """)
        st.sidebar.write("""
                ## DNA Tools
                """)
        sequence_file_path = st.sidebar.text_input(
            "DNA Input",
            "data_samples/Covid_sequence-NC_045512.fasta")
        if st.sidebar.button('Open'):
            sequence_data = dna_data.DNA()
            sequence_data.parse_dna_record(sequence_file_path)
            self.analyze_dna_sequence(sequence_data)

        self.sidebar_chemical_compound_menu()

    def sidebar_chemical_compound_menu(self):
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


ui_instance = ui_instance()