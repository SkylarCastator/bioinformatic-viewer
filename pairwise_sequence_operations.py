from scipy.spatial.distance import hamming
from skbio import DNA, Protein
import numpy as np


def return_hamming_distance_between_sequences(reference_seq, query_seq):
    reference_protein = skbio.protein(reference_seq)
    query_protein = skbio.protein(query_seq)
    return hamming(reference_protein, query_protein)


def show_difference_between_sequences(seq_1, seq_2):
    num_rows = len(seq_2)
    num_cols = len(seq_1)
    data = np.zeros(shape=(num_rows, num_cols), dtype=np.int)
    for row_number, row_character in enumerate(seq_2):
        for col_number, col_character in enumerate(seq_1):
            if row_character == col_character:
                data[row_number, col_number] = 1


def identify_longest_diagnals(data):
    summed_data = data.copy()
    for i in range(1, summed_data.shape[0]):
        for j in range(1, summed_data.shape[1]):
            if summed_data[i, j] > 0:
                summed_data[i, j] += summed_data[i-1][j-1]
    longest_diagonal_length = summed_data.max()
    #show


def needleman_wunsch_alignment(seq_1, seq_2):
    seq1 = Protein(seq_1)
    seq2 = Protein(seq_2)
    num_rows = len(seq1) + 1
    num_cols = len(seq2) + 1
    f_data = np.zeros(shape=(num_rows, num_cols), dtype=np.int)
    t_data = np.full(shape=(num_rows, num_cols), fill_value=" ", dtype=np.str)
    #show F and show T
    gap_penalty = 8     # d
    f_data[0][0] = 0
    for i in range(1, num_rows):
        f_data[i][0] = f_data[i-1][0] - gap_penalty
    for j in range(1, num_cols):
        f_data[0][j] = f_data[0][j-1] - gap_penalty
    # show F

    t_data[0][0] = "•"
    for i in range(1, num_rows):
        t_data[i][0] = "↑"
    for j in range(1, num_cols):
        t_data[0][j] = "←"
