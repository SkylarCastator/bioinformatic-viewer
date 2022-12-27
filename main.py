import streamlit as st
import sequence_analyzer as covid


class ui_instance:
    def __init__(self):
        st.write("""
        #Bio-Informatic reader
        Testing a user interface to view information of medicine and diseases.
        """)
        self.tools_sidebar()

    def tools_sidebar(self):
        add_selectbox = st.sidebar.selectbox(
            'How would you like to be contacted?',
            ('Email', 'Home phone', 'Mobile phone')
        )
        add_slider = st.sidebar.slider(
            'Select a range of values',
            0.0, 100.0, (25.0, 75.0)
        )


#Id NC_045512.2
#https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=graph
#https://www.rcsb.org/structure/6ZCO
#nv.demo()
#https://william-dawson.github.io/using-py3dmol.html
#https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753
sequence_file_path = "../data_samples/Covid_sequence-NC_045512.fasta"
sequence_data = covid.SequenceAnalyzer(sequence_file_path)
ui_instance = ui_instance()