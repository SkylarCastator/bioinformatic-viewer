import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

sequence_string = 'ATGTCTCGTAA'
my_seq = Seq(sequence_string)
print(my_seq)
print(my_seq.alphabet)

iupac_seq = Seq(sequence_string, IUPAC.unambiguous_dna)


print(len(my_seq))
print(my_seq[0:6])
print(my_seq.index('T'))

import matplotlib.pyplot as plt
from collections import Counter
dna_freq = Counter(my_seq)

#plt.bar(dna_freq.keys(), dna_freq.values())
#plt.show()


#Protien Synthesis
#A T -  2 hydrogen bonds
#G C - 3 hydorgen bonds
my_seq_complement = my_seq.complement()
print(my_seq_complement)
my_reverse_seq_complement = my_seq.reverse_complement()
print(my_reverse_seq_complement)

#transcription (DNA->RNA)
mRNA = my_seq.transcribe()
#translation to protiens (mRNA->Protein/Amino Acid)
proteins = mRNA.translate()

#Method 2 DNA->Amino Acids
proteins_2 = my_seq.translate()

#print(proteins.back_transcribe() == my_seq)

from Bio.SeqUtils import *
print(seq3(proteins))

from Bio.Data import CodonTable
print(CodonTable.unambiguous_dna_by_name['Standard'])