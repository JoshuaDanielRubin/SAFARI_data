import random
from typing import List, Dict
from collections import defaultdict

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

def apply_deamination_mutations(reads: List[str], mutation_rate: float) -> List[str]:
    """
    Apply deamination mutations (C->T or G->A) to the reads at a specified rate.
    """
    mutated_reads = []
    for read in reads:
        mutated_read = []
        for base in read:
            if base == "C" and random.random() < mutation_rate:
                mutated_read.append("T")
            elif base == "G" and random.random() < mutation_rate:
                mutated_read.append("A")
            else:
                mutated_read.append(base)
        mutated_reads.append(''.join(mutated_read))
    return mutated_reads

# Function for mismatch analysis

def find_mismatched_kmers(reads: List[str], k: int, w: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]]) -> List[int]:
    """
    Identify k-mers that match the k-mer index but not the rymer index, and return the number of mismatches observed.
    """
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
    return mismatch_counts

# Main Execution

sequence = read_fasta("rCRS.fa")
k = 8
w = 10

# Create minimizer and rymer tables
minimizer_table = create_index_table(sequence, k, w)
rymer_table = create_index_table(rymer_transform(sequence), k, w)

# Generate circular reads and introduce mutations
read_length = 50
mutation_rate = 0.5
circular_reads = generate_circular_reads(sequence, read_length, 5)
mutated_reads = apply_deamination_mutations(circular_reads, mutation_rate)

# Extract mismatch counts for the subset of k-mers
mismatch_counts = find_mismatched_kmers(mutated_reads, k, w, minimizer_table, rymer_table)

# Output can be plotted or further analyzed
print(mismatch_counts)
