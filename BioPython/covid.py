from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight, MeltingTemp as mt
from Bio.SeqUtils import GC, GC123,GC_skew
import pandas as pd
import py3Dmol
import nglview as nv
from Bio.PDB import PDBParser, MMCIFParser


class SequenceAnalyzer:
    def __init__(self, data_path):
        sequence_record = self.parse_record(data_path)
        dna, rna, proteins = self.retrieve_dna_rna_proteins(sequence_record)
        print(len(dna))
        #self.show_dna_count_plot(dna)
        print(molecular_weight(dna))
        print(GC(dna))
        print(self.at_content_percentage(dna))
        #Leading or lagging strand
        print(GC_skew(dna))
        print(mt.Tm_GC(dna, strict=False))
        #self.show_dna_count_plot(proteins)
        #Longest sequence
        proteins_split = proteins.split('*')
        df = pd.DataFrame({'amino_acids':proteins_split})
        df['count'] = df['amino_acids'].str.len()
        df.head()
        df.nlargest(10, 'count')
        self.show_3d_diagram()


    def retrieve_dna_rna_proteins(self, record_data):
        dna = record_data.seq
        print(dna)
        rna = dna.transcribe()
        print(rna)
        proteins = rna.translate()
        print(proteins)
        return dna, rna, proteins

    def parse_record(self, path):
        record_content = SeqIO.read(path, "fasta")
        # print(record_content)
        return record_content

    def at_content_percentage(self, dna):
        return (dna.count('A') + dna.count('T')/len(dna))*100

    def show_dna_count_plot(self, dna):
        n_count = Counter(dna)
        plt.bar(n_count.keys(), n_count.values())
        plt.show()

    def parse_3d_structure(self, pdb_path):
        parser = PDBParser()
        structure = parser.get_structure('PDB_structure_name', pdb_path)
        return structure

    def log_structure(self, structure):
        # strucutre -> model -> chain -> residue -> atom
        model = structure[0]
        for chain in model:
            print(f"Chain : {chain} as {chain.id}")
            for residue in chain:
                print(f"residue : {residue}")
                for atom in residue:
                    print(f"atom : {atom}")

    def show_3d_diagram(self):
        with open("6zco.pdb") as ifile:
            system = "".join([x for x in ifile])
        view = py3Dmol.view(width=400, height=300)
        view.addModelsAsFrames(system)
        view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
        view.zoomTo()
        view.show()
        #view = nv.show_biopython(structure)


#Id NC_045512.2
#https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=graph
#https://www.rcsb.org/structure/6ZCO
#nv.demo()
#https://william-dawson.github.io/using-py3dmol.html
#https://towardsdatascience.com/molecular-visualization-in-streamlit-using-rdkit-and-py3dmol-part-2-657d28152753
sequence_file_path = "Covid_sequence-NC_045512.fasta"
sequence_data = SequenceAnalyzer(sequence_file_path)





