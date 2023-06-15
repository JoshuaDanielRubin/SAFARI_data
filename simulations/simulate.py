import random
from collections import defaultdict
import argparse

def generate_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def reverse_complement(kmer):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(kmer))

def create_minimizer_index(S, W, k):
    minimizer_index = defaultdict(list)
    rymer_index = defaultdict(list)

    for i in range(len(S) - W + 1):
        window = S[i:i+W]
        for j in range(len(window) - k + 1):
            kmer = window[j:j+k]
            minimizer = min(kmer, reverse_complement(kmer))
            rymer = ''.join(['R' if base in 'GA' else 'Y' for base in kmer])
            minimizer_index[minimizer].append((i + j, kmer))
            rymer_index[rymer].append((i + j, kmer))

    return minimizer_index, rymer_index

def generate_reads(S, N, L, delta):
    reads = []
    for _ in range(N):
        start = random.randint(0, len(S) - L)
        read = list(S[start:start+L])
        for i in range(len(read)):
            if random.random() < delta:
                read[i] = random.choice('ACGT')
        reads.append((''.join(read), start))
    return reads

def find_seeds(index, read, k):
    seeds = defaultdict(list)
    total_seeds = 0
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        if kmer in index:
            seeds[kmer].extend(index[kmer])
            total_seeds += len(index[kmer])
    return seeds, total_seeds

def check_seeds(seeds, read_origin):
    correct = 0
    for kmer, locations in seeds.items():
        for location in locations:
            if location[0] == read_origin:
                correct += 1
    return correct

def main(N, delta, k, W):
    S = generate_dna(10000)
    minimizer_index, rymer_index = create_minimizer_index(S, W, k)
    L = 50
    reads = generate_reads(S, N, L, delta)

    total_correct_minimizer_seeds = 0
    total_correct_rymer_seeds = 0
    total_minimizer_seeds = 0
    total_rymer_seeds = 0

    for read, origin in reads:
        rymer_read = ''.join(['R' if base in 'GA' else 'Y' for base in read])
        minimizer_seeds, minimizer_seeds_count = find_seeds(minimizer_index, read, k)
        rymer_seeds, rymer_seeds_count = find_seeds(rymer_index, rymer_read, k)

        correct_minimizer_seeds = check_seeds(minimizer_seeds, origin)
        correct_rymer_seeds = check_seeds(rymer_seeds, origin)
        total_correct_minimizer_seeds += correct_minimizer_seeds
        total_correct_rymer_seeds += correct_rymer_seeds
        total_minimizer_seeds += minimizer_seeds_count
        total_rymer_seeds += rymer_seeds_count

    fraction_correct_minimizer_seeds = total_correct_minimizer_seeds / total_minimizer_seeds
    fraction_correct_rymer_seeds = total_correct_rymer_seeds / total_rymer_seeds

    print(f'Total correct minimizer seeds: {total_correct_minimizer_seeds}  Sensitivity of minimizer seeds: {fraction_correct_minimizer_seeds:.4f}  Precision of minimizer seeds: {total_correct_minimizer_seeds / total_minimizer_seeds:.4f}')
    print(f'Total correct rymer seeds: {total_correct_rymer_seeds}  Sensitivity of rymer seeds: {fraction_correct_rymer_seeds:.4f}  Precision of rymer seeds: {total_correct_rymer_seeds / total_rymer_seeds:.4f}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=100, help='Number of reads')
    parser.add_argument('--delta', type=float, default=0.01, help='Mutation rate')
    parser.add_argument('--k', type=int, default=5, help='Length of k-mer')
    parser.add_argument('--W', type=int, default=10, help='Length of window')
    args = parser.parse_args()

    assert 0 <= args.delta <= 1, "Mutation rate must be between 0 and 1"
    assert args.N > 0, "Number of reads must be greater than 0"
    assert args.k > 0, "Length of k-mer must be greater than 0"
    assert args.W >= args.k, "Length of window must be greater than or equal to length of k-mer"

    main(args.N, args.delta, args.k, args.W)

