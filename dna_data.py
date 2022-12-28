from Bio import SeqIO
from Bio.SeqUtils import molecular_weight, MeltingTemp as mt
from Bio.SeqUtils import GC
import neatbio.sequtils as utils


class DNA:
    def __init__(self):
        self.name = ""
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
        self.name = ""
        self.dna = sequence_record.seq
        self.rna = self.dna.transcribe()
        self.proteins = self.rna.translate()
        proteins_three_letter = utils.convert_lto3(str(self.proteins).replace("*", ""))
        proteins_full_names = utils.get_acid_name(proteins_three_letter)
        self.dna_length = len(self.dna)
        self.molecular_weight = molecular_weight(self.dna)
        self.a_t_percentage = self.a_t_content_percentage(self.dna)
        self.g_c_percentage = GC(self.dna)
        self.melting_point_gc = str(mt.Tm_GC(self.dna, strict=False)) + "C"
        self.melting_point_wallace = str(mt.Tm_Wallace(self.dna)) + "C"

    def a_t_content_percentage(self, dna):
        return (dna.count('A') + dna.count('T')/len(dna))*100

    def parse_record(self, path):
        record_content = SeqIO.read(path, "fasta")
        return record_content