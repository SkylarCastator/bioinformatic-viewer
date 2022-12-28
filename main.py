import streamlit as st
import py3Dmol
from stmol import showmol
import pandas as pd
import numpy as np
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

    def show_3d_mol(self, protein, bcolor, spin, style):
        xyzview = py3Dmol.view(query='pdb:' + protein)
        xyzview.setStyle({style: {'color': 'spectrum'}})
        xyzview.setBackgroundColor(bcolor)
        if spin:
            xyzview.spin(True)
        else:
            xyzview.spin(False)
        xyzview.zoomTo()
        showmol(xyzview, height=500, width=800)

    def analyze_dna_sequence(self, dna_sequence):
        self.show_3d_mol("1A2C", '#ffffff', True, 'cartoon')
        d = {
            'Key': [
                'Name',
                'Length',
                'Molecular Weight',
                'GC Percentage',
                'AT Percentage',
                'Melting Point GC',
                'Melting Point Wallace'],
            'Value': [
                dna_sequence.name,
                dna_sequence.dna_length,
                dna_sequence.molecular_weight,
                dna_sequence.g_c_percentage,
                dna_sequence.a_t_percentage,
                dna_sequence.melting_point_gc,
                dna_sequence.melting_point_wallace]
        }
        df = pd.DataFrame(data=d)
        st.dataframe(df)

        st.subheader('Nitrogenous Bases')
        n_count = Counter(dna_sequence.dna)
        df = pd.DataFrame.from_dict(n_count, orient='index')
        st.bar_chart(df)

        st.subheader('Proteins')
        n_count = Counter(dna_sequence.proteins)
        df = pd.DataFrame.from_dict(n_count, orient='index')
        st.bar_chart(df)

    def show_table_view(self):
        df = pd.DataFrame(
            np.random.randn(10, 5),
            colums=(' col %d' % i for i in range(5))
        )
        st.table(df)


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