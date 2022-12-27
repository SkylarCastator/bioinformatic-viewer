from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight, MeltingTemp as mt
from Bio.SeqUtils import GC, GC123, GC_skew
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser


class DNA:
    def __init__(self):
        self.dna_length = ""
        self.molecular_weight = ""
        self.g_c_percentage = ""
        self.a_t_percentage = ""
        self.melting_point_gc = ""
        self.melting_point_wallace = ""
        self.dna = None
        self.rna = None
        self.proteins = None

    def parse_dna_record(self, data_path):
        sequence_record = self.parse_record(data_path)
        self.dna = sequence_record.seq
        self.rna = self.dna.transcribe()
        self.proteins = self.rna.translate()
        self.dna_length = len(self.dna)
        self.molecular_weight = molecular_weight(self.dna)
        self.a_t_percentage = self.a_t_content_percentage(self.dna)
        self.g_c_percentage = GC(self.dna)
        self.melting_point_gc = str(mt.Tm_GC(self.dna, strict=False)) + "C"
        self.melting_point_wallace = str(mt.Tm_Wallace(self.dna)) + "C"


        # Leading or lagging strand
        #print(GC_skew(dna))

        # Longest sequence
        #proteins_split = proteins.split('*')
        #df = pd.DataFrame({'amino_acids': proteins_split})
        #df['count'] = df['amino_acids'].str.len()
        #df.head()
        #df.nlargest(10, 'count')

    def a_t_content_percentage(self, dna):
        return (dna.count('A') + dna.count('T')/len(dna))*100

    def parse_record(self, path):
        record_content = SeqIO.read(path, "fasta")
        return record_content