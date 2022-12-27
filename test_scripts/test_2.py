from Bio.SeqUtils import GC, MeltingTemp as mt
from Bio.Seq import Seq

seq = 'GCTGTCGTCATT'
print(GC(seq))
print(str(mt.Tm_Wallace(seq)) + "C")
print(str(mt.Tm_GC(seq)) + "C")

# GC Skew
# indicates leading or lagging strain
# pos -> leading
# neg -> Lagging
from Bio.SeqUtils import GC123, GC_skew, nt_search

print(GC123(seq))
print(GC_skew(seq))

# Tinker Graphing:
# xGC_skew

# search subsequence
sub_seq = 'AAT'
nt_search(str(seq), str(sub_seq))

# Sequence Alignment
# Homology
# match, mismatch, gap
# Alignment type Local(Water) verse global(Needle)
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
seq1 = Seq('ACTCGT')
seq2 = Seq('ATTCG')
alignments = pairwise2.align.globalxx(seq1, seq2)
print(alignments)
print(format_alignment(*alignments[0]))
alignments = pairwise2.align.localxx(seq1, seq2, one_alignment_only=True, score_only=True)
print(alignments)
#print(format_alignment(*alignments[0]))

alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
print(alignments)