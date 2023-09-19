import random
from typing import List, Dict
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import gzip
from concurrent.futures import ProcessPoolExecutor
import os

# Create complement mapping outside of the function
COMPLEMENT = str.maketrans('ACTGNRY', 'TGACNYR')

# Add this function to read sequences from a gzipped FASTQ file
def read_fastq(file_path: str) -> List[str]:
    with gzip.open(file_path, 'rt') as f:
        sequences = []
        while True:
            f.readline()  # skip the name line
            seq = f.readline().strip()  # get the sequence line
            if not seq:
                break
            f.readline()  # skip the plus line
            f.readline()  # skip the quality line
            sequences.append(seq.upper())
    return sequences

# Functions for reading and preprocessing
def read_fasta(file_path: str) -> str:
    with open(file_path, 'r') as f:
        sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
    return sequence.upper()

def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

RYMER_MAP = {'C': 'C', 'T': 'C', 'G': 'A', 'A': 'A', 'N': 'N', 'R': 'R', 'Y': 'Y'}

def rymer_transform(seq: str) -> str:
    return ''.join(RYMER_MAP.get(base, 'N') for base in seq)

def create_minimizer(seq: str, k: int, w: int) -> str:
    return min(seq[i:i+k] for i in range(w - k + 1))

def create_index_table(sequence: str, k: int, w: int) -> Dict[str, List[int]]:
    table = defaultdict(list)
    for i in range(len(sequence) - w + 1):
        window = sequence[i:i+w]
        minimizer = create_minimizer(window, k, w)
        table[minimizer].append(i)
    return table

def process_kmer(args):
    kmer, i, rymer_set, sequence, k = args
    rymer = rymer_transform(kmer)
    rc_kmer = reverse_complement(kmer)
    rc_rymer = rymer_transform(rc_kmer)

    rymer_found = rymer in rymer_set or rc_rymer in rymer_set

    if rymer_found:
        ref_segment = sequence[i:i+k]

        mismatch_count = sum((a == 'T' and b == 'C') or (a == 'A' and b == 'G') for a, b in zip(kmer, ref_segment))
        exact_match = kmer == ref_segment
        return mismatch_count, exact_match
    else:
        return None, None

def process_read(args):
    read, k, rymer_set, sequence = args
    mismatch_counts = []
    total_kmers = 0
    exact_matches = 0

    kmer_args = [(read[i:i+k], i, rymer_set, sequence, k) for i in range(len(read) - k + 1)]

    with ProcessPoolExecutor(max_workers=60) as executor:
        results = executor.map(process_kmer, kmer_args)

    for res in results:
        mismatch_count, exact_match = res
        if mismatch_count is not None:
            mismatch_counts.append(mismatch_count)
            total_kmers += 1
            exact_matches += exact_match

    return mismatch_counts, exact_matches, total_kmers

def find_deamination_mismatches(read: str, k: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]], sequence: str) -> List[int]:
    rymer_set = set(rymer_table.keys())
    args = (read, k, rymer_set, sequence)
    mismatch_counts, exact_matches, total_kmers = process_read(args)
    return mismatch_counts, exact_matches, total_kmers

# Main code for generating the plot
sequence = read_fasta("rCRS.fa")
bacterial_reference = read_fasta("refSoil.fa")
#bacterial_reference = read_fasta("small.fa")
print("read the references")
k_values = list(range(3, 16))
average_mismatches = []
exact_match_fractions = []

for k in k_values:
    print("on k = " + str(k))
    w = k + 2
    minimizer_table = create_index_table(sequence, k, w)
    print("minimizer index created")
    rymer_table = create_index_table(rymer_transform(sequence), k, w)
    print("rymer index created")
    mismatch_counts, exact_matches, total_kmers = find_deamination_mismatches([bacterial_reference], k, minimizer_table, rymer_table, sequence)
    print("deam mismatches found")
    non_zero_mismatches = [count for count in mismatch_counts if count > 0]
    average_mismatch = sum(non_zero_mismatches) / len(non_zero_mismatches) if non_zero_mismatches else 0
    average_mismatches.append(average_mismatch / k)
    exact_match_fractions.append(exact_matches / total_kmers if total_kmers else 0)

plt.figure(figsize=(10, 5))

# First subplot
plt.subplot(1, 2, 1)
plt.plot(k_values, average_mismatches, marker='o', linestyle='-')
plt.xticks(k_values)
plt.xlabel('Value of k')
plt.ylabel('Mismatch Proportion')
plt.title('Sequence Similarity as a Function of k (Deamination)')
plt.grid(True)

# Curve fitting for first subplot
def power_law(x, a, b):
    return a * np.power(x, b)

params, _ = curve_fit(power_law, k_values, average_mismatches)
a, b = params

# Annotate the plot with the fitted parameters
annotation_text = f'a={a:.4f}, b={b:.4f}'

x_fit = np.linspace(min(k_values), max(k_values), 1000)
y_fit = power_law(x_fit, *params)
plt.plot(x_fit, y_fit, label='Power-law fit', linestyle='--')
plt.legend()

# Second subplot
plt.subplot(1, 2, 2)
plt.plot(k_values, exact_match_fractions, marker='o', linestyle='-', color='red')
plt.xticks(k_values)
plt.xlabel('Value of k')
plt.ylabel('Exact Match Fraction')
plt.title('Exact Match Fraction as a Function of k (Deamination)')
plt.grid(True)

plt.subplot(1, 2, 1)  # Make sure we're on the first subplot
plt.annotate(f'a={a:.4f}, b={b:.4f}', xy=(0.6, 0.2), xycoords='axes fraction')
plt.tight_layout()
plt.savefig("mismatch.png")

# Print the parameters
print(f"Fitted parameters: a = {a}, b = {b}")

