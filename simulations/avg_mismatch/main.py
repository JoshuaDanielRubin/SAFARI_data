import random
from typing import List, Dict
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Functions for reading and preprocessing
def read_fasta(file_path: str) -> str:
    with open(file_path, 'r') as f:
        sequence = ''.join(line.strip() for line in f.readlines()[1:])
    return sequence.upper()

def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

# Functions for index creation
def rymer_transform(seq: str) -> str:
    return ''.join(['C' if base in ['C', 'T'] else 'A' if base in ['G', 'A'] else 'N' for base in seq])

def create_minimizer(seq: str, k: int, w: int) -> str:
    return min(seq[i:i+k] for i in range(w - k + 1))

def create_index_table(sequence: str, k: int, w: int) -> Dict[str, List[int]]:
    table = defaultdict(list)
    for i in range(len(sequence) - w + 1):
        window = sequence[i:i+w]
        minimizer = create_minimizer(window, k, w)
        table[minimizer].append(i)
    return table

# Functions for read generation and mutation introduction
def fragment_genome(sequence: str, mean_fragment_size: int, std_dev: int, num_fragments: int) -> List[str]:
    fragments = []
    for _ in range(num_fragments):
        start = random.randint(0, len(sequence) - mean_fragment_size)
        fragment_size = int(np.random.normal(mean_fragment_size, std_dev))
        fragment = sequence[start:start+fragment_size]
        fragments.append(fragment)
    return fragments

def generate_circular_reads(sequence: str, read_length: int, num_reads: int = None) -> List[str]:
    assert len(sequence) >= read_length, "Sequence length should be greater than or equal to read_length."
    circular_sequence = sequence + sequence[:read_length-1]
    if num_reads is None:
        num_reads = len(sequence)
    reads = [circular_sequence[i:i+read_length] for i in range(num_reads)]
    assert all(len(read) == read_length for read in reads), "All reads should have the same length."
    return reads

def briggs_model_damage_probability(x, P0=0.9, lambda_val=10):
    import math
    return P0 * math.exp(-x/lambda_val)

def apply_deamination_mutations(reads: List[str], mutation_rate: float) -> List[str]:
    mutated_reads = []
    for read in reads:
        mutated_read = []
        for i, base in enumerate(read):
            prob_start = briggs_model_damage_probability(i)
            prob_end = briggs_model_damage_probability(len(read) - 1 - i)
            prob = max(prob_start, prob_end)
            if base == 'C' and random.random() < prob * mutation_rate:
                mutated_read.append('T')
            else:
                mutated_read.append(base)
        mutated_reads.append("".join(mutated_read))
    return mutated_reads

# Modified function to count only deamination-specific mismatches
def find_deamination_mismatches(reads: List[str], k: int, w: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]], sequence: str) -> List[int]:
    total_kmers = 0
    exact_matches = 0
    mismatch_counts = []
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            rymer = rymer_transform(kmer)
            rc_kmer = reverse_complement(kmer)
            rc_rymer = rymer_transform(rc_kmer)
            kmer_found = kmer in minimizer_table or rc_kmer in minimizer_table
            rymer_found = rymer in rymer_table or rc_rymer in rymer_table
            if rymer_found:
                ref_segment = sequence[i:i+k]
                mismatch_count = sum(
                    1 for a, b in zip(kmer, ref_segment) if (a == 'T' and b == 'C') or (a == 'A' and b == 'G')
                )
                mismatch_counts.append(mismatch_count)
                total_kmers += 1
                if kmer == ref_segment:
                    exact_matches += 1
    return mismatch_counts, exact_matches, total_kmers

# Main code for generating the plot
sequence = read_fasta("rCRS.fa")
k_values = list(range(4, 15))
average_mismatches = []
exact_match_fractions = []

for k in k_values:
    w = k + 2
    minimizer_table = create_index_table(sequence, k, w)
    rymer_table = create_index_table(rymer_transform(sequence), k, w)
    fragments = fragment_genome(sequence, 150, 0, 200)
    all_reads = []
    for fragment in fragments:
        reads_from_fragment = generate_circular_reads(fragment, 75, 1)
        all_reads.extend(reads_from_fragment)
    mutated_reads = apply_deamination_mutations(all_reads, 0.2)
    mismatch_counts, exact_matches, total_kmers = find_deamination_mismatches(mutated_reads, k, w, minimizer_table, rymer_table, sequence)
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
plt.ylabel('Sequence Similarity')
plt.title('Sequence Similarity as a Function of k (Deamination)')
plt.grid(True)

# Curve fitting for first subplot
def power_law(x, a, b):
    return a * np.power(x, b)

params, _ = curve_fit(power_law, k_values, average_mismatches)
a, b = params

# Annotate the plot with the fitted parameters
annotation_text = f'a={a:.4f}, b={b:.4f}'
#plt.annotate(annotation_text, xy=(0.6, 0.2), xycoords='axes fraction')

x_fit = np.linspace(min(k_values), max(k_values), 500)
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

plt.tight_layout()
plt.savefig("mismatch.png")

# Print the parameters
print(f"Fitted parameters: a = {a}, b = {b}")
