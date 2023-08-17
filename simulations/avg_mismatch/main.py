import random
from typing import List, Dict
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

# Functions for reading and preprocessing

def read_fasta(file_path: str) -> str:
    """
    Read a FASTA file and return the sequence as a single string.
    """
    with open(file_path, 'r') as f:
        sequence = ''.join(line.strip() for line in f.readlines()[1:])
    return sequence.upper()

def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence, safely handling 'N' bases.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

# Functions for index creation

def rymer_transform(seq: str) -> str:
    """
    Transforms a DNA sequence into its rymer equivalent.
    C and T become C.
    G and A become A.
    """
    return ''.join(['C' if base in ['C', 'T'] else 'A' if base in ['G', 'A'] else 'N' for base in seq])

def create_minimizer(seq: str, k: int, w: int) -> str:
    """
    Extract the lexicographically smallest k-mer from a window of length w in the sequence.
    """
    return min(seq[i:i+k] for i in range(w - k + 1))

def create_index_table(sequence: str, k: int, w: int) -> Dict[str, List[int]]:
    """
    Create a minimizer index table for a sequence.
    """
    table = defaultdict(list)
    for i in range(len(sequence) - w + 1):
        window = sequence[i:i+w]
        minimizer = create_minimizer(window, k, w)
        table[minimizer].append(i)
    return table

# Functions for read generation and mutation introduction

def fragment_genome(sequence: str, mean_fragment_size: int, std_dev: int, num_fragments: int) -> List[str]:
    """
    Fragment the genome into smaller pieces.
    The size of the fragments is drawn from a normal distribution with the specified mean and standard deviation.
    """
    fragments = []
    for _ in range(num_fragments):
        start = random.randint(0, len(sequence) - mean_fragment_size)
        fragment_size = int(np.random.normal(mean_fragment_size, std_dev))
        fragment = sequence[start:start+fragment_size]
        fragments.append(fragment)
    return fragments

def generate_circular_reads(sequence: str, read_length: int, num_reads: int = None) -> List[str]:
    """
    Generate circular reads of the specified length from the given sequence.
    """
    # Ensure that the sequence length is at least equal to the read_length
    assert len(sequence) >= read_length, "Sequence length should be greater than or equal to read_length."

    # Create a circular sequence by appending the start of the sequence to its end
    circular_sequence = sequence + sequence[:read_length-1]

    # If num_reads is not provided, default to the length of the sequence
    if num_reads is None:
        num_reads = len(sequence)

    # Generate the reads
    reads = [circular_sequence[i:i+read_length] for i in range(num_reads)]

    # Assert all reads are of the same length
    assert all(len(read) == read_length for read in reads), "All reads should have the same length."

    return reads

# Implementing the Briggs model for ancient DNA damage
def briggs_model_damage_probability(x, P0=0.9, lambda_val=10):
    """Calculate the probability of damage at position x using the Briggs model."""
    import math
    return P0 * math.exp(-x/lambda_val)

# Incorporating the Briggs model into the mutation introduction function
def apply_deamination_mutations(reads: List[str], mutation_rate: float) -> List[str]:
    """
    Introduce deamination mutations to a list of reads using the Briggs model.
    C -> T mutations at the start and end of fragments.
    """
    mutated_reads = []
    
    for read in reads:
        mutated_read = []
        
        for i, base in enumerate(read):
            # Calculate probability of mutation at this position using the Briggs model
            # Consider both ends (start and end) of the fragment
            prob_start = briggs_model_damage_probability(i)
            prob_end = briggs_model_damage_probability(len(read) - 1 - i)
            prob = max(prob_start, prob_end)
            
            # Apply mutation with the calculated probability
            if base == 'C' and random.random() < prob * mutation_rate:
                mutated_read.append('T')
            else:
                mutated_read.append(base)
        
        mutated_reads.append("".join(mutated_read))
    
    return mutated_reads

# Function for mismatch analysis

def find_mismatched_kmers(reads: List[str], k: int, w: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]]) -> List[int]:
    """
    Identify k-mers that match the k-mer index but not the rymer index, and return the number of mismatches observed.
    """
    total_kmers=0
    exact_matches=0

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
                # Fix: Directly get the corresponding sequence segment instead of using a different position
                ref_segment = sequence[i:i+k]
                mismatch_count = sum(1 for a, b in zip(kmer, ref_segment) if (a == 'T' and b == 'C') or (a == 'A' and b == 'G'))
                mismatch_counts.append(mismatch_count)

                # Check for exact match and increment counters
                total_kmers += 1
                if kmer == ref_segment:
                    exact_matches += 1

    return mismatch_counts, exact_matches, total_kmers

def compute_mismatch_average_for_k(sequence, k):
    w = k + 2
    
    # Create minimizer and rymer tables
    minimizer_table = create_index_table(sequence, k, w)
    rymer_table = create_index_table(rymer_transform(sequence), k, w)

    # Generate circular reads and introduce mutations
    read_length = 75
    mutation_rate = 0.2
    mean_fragment_size = 150
    std_dev = 0
    num_fragments = 1000
    num_reads_per_fragment = 2

    fragments = fragment_genome(sequence, mean_fragment_size, std_dev, num_fragments)
    all_reads = []
    for fragment in fragments:
        reads_from_fragment = generate_circular_reads(fragment, read_length, num_reads_per_fragment)
        all_reads.extend(reads_from_fragment)

    mutated_reads = apply_deamination_mutations(all_reads, mutation_rate)
    mismatch_counts, exact_matches, total_kmers = find_mismatched_kmers(mutated_reads, k, w, minimizer_table, rymer_table)

    # Calculate the average excluding zero mismatches
    non_zero_mismatches = [count for count in mismatch_counts if count > 0]
    average_mismatch = sum(non_zero_mismatches) / len(non_zero_mismatches) if non_zero_mismatches else 0

    return average_mismatch, exact_matches / total_kmers if total_kmers else 0


sequence = read_fasta("rCRS.fa")
# Define range of k values
k_values = list(range(4, 31))  # Modify as per your needs
average_mismatches = []

# Looping through k values
exact_match_fractions = []
for k in k_values:
    print(k)
    avg_mismatch, exact_match_fraction = compute_mismatch_average_for_k(sequence, k)
    avg_mismatch_normalized = avg_mismatch / k
    average_mismatches.append(avg_mismatch_normalized)
    exact_match_fractions.append(exact_match_fraction)

# Plotting the results for average mismatches
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(k_values, average_mismatches, marker='o', linestyle='-')
plt.xticks(k_values)
plt.xlabel('Value of k')
plt.ylabel('Sequence Similarity')
plt.title('Sequence Similarity as a Function of k')
plt.grid(True)

# Plotting the results for exact match fractions
plt.subplot(1, 2, 2)
plt.plot(k_values, exact_match_fractions, marker='o', linestyle='-', color='red')
plt.xticks(k_values)
plt.xlabel('Value of k')
plt.ylabel('Exact Match Fraction')
plt.title('Exact Match Fraction as a Function of k')
plt.grid(True)

plt.tight_layout()
plt.savefig("mismatch.png")
plt.close()
