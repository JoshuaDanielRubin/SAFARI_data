import argparse
import random
from Bio import SeqIO
from collections import defaultdict
from sklearn import metrics
import numpy as np

# Function to convert a sequence to RY representation
def to_ry(sequence):
    conversion = {'A': 'R', 'G': 'R', 'C': 'Y', 'T': 'Y'}
    return ''.join(conversion[base] for base in sequence if base in 'ACGT')

# Function to extract k-mers from a sequence
def get_kmers(sequence, k):
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    valid_bases = set('ACGTRY')
    kmers = [kmer for kmer in kmers if set(kmer).issubset(valid_bases)]
    return kmers

# Function to construct an index from a list of k-mers
def construct_index(kmers):
    index = defaultdict(list)
    for i, kmer in enumerate(kmers):
        index[kmer].append(i)
    return index

# Function to simulate ancient DNA damage
def simulate_damage(sequence, rate):
    conversion = {'C': 'T', 'G': 'A', 'A': 'A', 'T': 'T'}
    damaged_sequence = ''
    for base in sequence:
        if base in 'CG' and random.random() < rate:
            damaged_sequence += conversion[base]
        else:
            damaged_sequence += base
    return damaged_sequence

# Function to simulate k-mer and RY-mer finding in indexes and compute sensitivity, specificity, and AUC
def simulate(index_kmer, index_ry, sequence, k):
    true_positive_kmer = 0
    false_positive_kmer = 0

    true_positive_ry = 0
    false_positive_ry = 0

    y_true_kmer = []
    y_pred_kmer = []

    y_true_ry = []
    y_pred_ry = []

    for start in range(len(sequence) - k):
        kmer = sequence[start:start + k]
        ry_kmer = to_ry(kmer)

        if kmer in index_kmer:
            if start in index_kmer[kmer]:
                true_positive_kmer += 1
                y_pred_kmer.append(1)
            else:
                false_positive_kmer += 1
                y_pred_kmer.append(0)
            y_true_kmer.append(1)
        else:
            y_pred_kmer.append(0)
            y_true_kmer.append(0)

        if ry_kmer in index_ry:
            if start in index_ry[ry_kmer]:
                true_positive_ry += 1
                y_pred_ry.append(1)
            else:
                false_positive_ry += 1
                y_pred_ry.append(0)
            y_true_ry.append(1)
        else:
            y_pred_ry.append(0)
            y_true_ry.append(0)

    sensitivity_kmer = true_positive_kmer / (len(sequence) - k)
    specificity_kmer = (len(sequence) - k - false_positive_kmer) / (len(sequence) - k)
    auc_kmer = metrics.roc_auc_score(y_true_kmer, y_pred_kmer)

    sensitivity_ry = true_positive_ry / (len(sequence) - k)
    specificity_ry = (len(sequence) - k - false_positive_ry) / (len(sequence) - k)
    auc_ry = metrics.roc_auc_score(y_true_ry, y_pred_ry)

    return (sensitivity_kmer, specificity_kmer, auc_kmer), (sensitivity_ry, specificity_ry, auc_ry)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Generate k-mer and RY-mer indexes from a genome')
    parser.add_argument('--genome', type=str, default='rCRS.fa', help='Path to the reference genome fasta file')
    parser.add_argument('--k', type=int, default=10, help='Size of the k-mers')
    parser.add_argument('--damage_rate', type=float, default=0.1, help='Rate of ancient DNA damage (C->T and G->A mutations)')
    args = parser.parse_args()

    # Load the reference genome
    genome = SeqIO.read(args.genome, 'fasta')
    genome_str = str(genome.seq)

    # Simulate ancient DNA damage
    damaged_genome_str = simulate_damage(genome_str, args.damage_rate)

    # Generate k-mers and construct k-mer index
    kmers = get_kmers(damaged_genome_str, args.k)
    kmer_index = construct_index(kmers)

    # Convert genome to RY representation, generate RY-mers and construct RY-mer index
    ry_genome_str = to_ry(damaged_genome_str)
    rymers = get_kmers(ry_genome_str, args.k)
    rymer_index = construct_index(rymers)

    # Run simulation for k-mer and RY-mer indexes and compute sensitivity, specificity, and AUC
    (sensitivity_kmer, specificity_kmer, auc_kmer), (sensitivity_ry, specificity_ry, auc_ry) = simulate(kmer_index, rymer_index, damaged_genome_str, args.k)

    print(f'For k-mer index: Sensitivity = {sensitivity_kmer}, Specificity = {specificity_kmer}, AUC = {auc_kmer}')
    print(f'For RY-mer index: Sensitivity = {sensitivity_ry}, Specificity = {specificity_ry}, AUC = {auc_ry}')

if __name__ == "__main__":
    main()

