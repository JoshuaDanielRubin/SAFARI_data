import argparse
import random
from Bio import SeqIO
from collections import defaultdict

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

def simulate(index, num_kmers, k, damage_rate):
    found = 0
    for _ in range(num_kmers):
        kmer = ''.join(random.choice('ACGT') for _ in range(k))
        if random.random() < damage_rate:
            kmer = simulate_damage(kmer, damage_rate)
        if kmer in index:
            found += 1
    return found

# Function to simulate RY-mer finding in an index
def simulate_ry(index, num_kmers, k, damage_rate):
    found = 0
    for _ in range(num_kmers):
        kmer = ''.join(random.choice('ACGT') for _ in range(k))
        kmer = to_ry(kmer)
        if kmer in index:
            found += 1
    return found

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

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Generate k-mer and RY-mer indexes from a genome')
    parser.add_argument('--genome', type=str, default='rCRS.fa', help='Path to the reference genome fasta file')
    parser.add_argument('--k', type=int, default=10, help='Size of the k-mers')
    parser.add_argument('--n', type=int, default=10000, help='Number of k-mers for simulation')
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

    # Run simulation for k-mer index
    found_kmers = simulate(kmer_index, args.n, args.k, args.damage_rate)

    print(f'Found {found_kmers} out of {args.n} random k-mers in k-mer index')

    # Run simulation for RY-mer index
    found_rymers = simulate_ry(rymer_index, args.n, args.k, args.damage_rate)
    print(f'Found {found_rymers} out of {args.n} random k-mers in RY-mer index')

if __name__ == "__main__":
    main()

