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
    read, k, w, minimizer_set, rymer_set, sequence = args
    mismatch_counts = []
    total_kmers = 0
    exact_matches = 0

    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        rymer = rymer_transform(kmer)
        rc_kmer = reverse_complement(kmer)
        rc_rymer = rymer_transform(rc_rymer)

        rymer_found = rymer in rymer_set or rc_rymer in rymer_set

        if rymer_found:
            ref_segment = sequence[i:i+k]

            mismatch_count = sum((a == 'T' and b == 'C') or (a == 'A' and b == 'G') for a, b in zip(kmer, ref_segment))
            mismatch_counts.append(mismatch_count)
            total_kmers += 1

            if kmer == ref_segment:
                exact_matches += 1

    return mismatch_counts, exact_matches, total_kmers

def find_deamination_mismatches(reads: List[str], k: int, w: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]], sequence: str) -> List[int]:
    total_kmers = 0
    exact_matches = 0
    mismatch_counts = []

    minimizer_set = set(minimizer_table.keys())
    rymer_set = set(rymer_table.keys())

    kmer_args = []

    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            kmer_args.append((kmer, i, rymer_set, sequence, k))

    # Subsample 0.1% of the kmers
    subsample_size = max(1, int(0.01 * len(kmer_args)))
    kmer_args_subsample = random.sample(kmer_args, subsample_size)

    for args in kmer_args_subsample:
        mc, em = process_kmer(args)
        if mc is not None and em is not None:
            mismatch_counts.append(mc)
            total_kmers += 1
            if em:
                exact_matches += 1

    return mismatch_counts, exact_matches, total_kmers

def process_k(k):
    print("k= " +str(k))
    w = k + 2
    minimizer_table = create_index_table(sequence, k, w)
    rymer_table = create_index_table(rymer_transform(sequence), k, w)
    mismatch_counts, exact_matches, total_kmers = find_deamination_mismatches([bacterial_reference], k, w, minimizer_table, rymer_table, sequence)
    non_zero_mismatches = [count for count in mismatch_counts if count > 0]
    average_mismatch = sum(non_zero_mismatches) / len(non_zero_mismatches) if non_zero_mismatches else 0
    average_mismatch /= k
    exact_match_fraction = exact_matches / total_kmers if total_kmers else 0
    return average_mismatch, exact_match_fraction

# Main code for generating the plot
sequence = read_fasta("rCRS.fa")
bacterial_reference = read_fasta("refSoilSmall.fa")

k_values = list(range(3, 31))
results = []

with ProcessPoolExecutor(max_workers=1) as executor:
    results = list(executor.map(process_k, k_values))

average_mismatches, exact_match_fractions = zip(*results)

# Plotting the results
plt.figure(figsize=(10, 5))

# Plotting the mismatch proportions
plt.plot(k_values, average_mismatches, marker='o', linestyle='-')
plt.xticks(k_values)
plt.xlabel('Value of k')
plt.ylabel('Mismatch Proportion')
plt.title('Sequence Similarity as a Function of k')
plt.grid(True)

# Curve fitting for the plot
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

plt.annotate(f'a={a:.4f}, b={b:.4f}', xy=(0.6, 0.2), xycoords='axes fraction')
plt.tight_layout()
plt.savefig("mismatch.png")

# Print the parameters
print(f"Fitted parameters: a = {a}, b = {b}")
