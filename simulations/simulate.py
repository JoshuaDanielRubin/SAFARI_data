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
            minimizer_index[minimizer].append(i + j)
            rymer_index[rymer].append(i + j)

    return minimizer_index, rymer_index

def generate_reads(S, N, L, delta):
    reads = []
    for _ in range(N):
        start = random.randint(0, len(S) - L)
        read = list(S[start:start+L])
        for i in range(len(read)):
            if random.random() < delta:
                if read[i] == 'C':
                    read[i] = 'T'
                elif read[i] == 'G':
                    read[i] = 'A'
        reads.append((''.join(read), start))
    return reads

def find_seeds(index, read, k, origin):
    seeds = []
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        if kmer in index:
            seeds.extend([(pos, pos == origin + i) for pos in index[kmer]])
    return seeds

def compute_stats(seeds):
    total_seeds_found = len(seeds)
    correct_seeds = sum(correct for _, correct in seeds)
    if total_seeds_found > 0:
        precision = correct_seeds / total_seeds_found
    else:
        precision = 0
    return total_seeds_found, precision

def main(N, L, delta, k, W):
    S = generate_dna(10000)
    minimizer_index, rymer_index = create_minimizer_index(S, W, k)
    reads = generate_reads(S, N, L, delta)

    total_minimizer_seeds = 0
    total_rymer_seeds = 0
    minimizer_precision = 0
    rymer_precision = 0

    for read, origin in reads:
        minimizer_seeds = find_seeds(minimizer_index, read, k, origin)
        rymer_read = ''.join(['R' if base in 'GA' else 'Y' for base in read])
        rymer_seeds = find_seeds(rymer_index, rymer_read, k, origin)

        minimizer_seeds_found, minimizer_seed_precision = compute_stats(minimizer_seeds)
        rymer_seeds_found, rymer_seed_precision = compute_stats(rymer_seeds)

        total_minimizer_seeds += minimizer_seeds_found
        total_rymer_seeds += rymer_seeds_found
        minimizer_precision += minimizer_seed_precision
        rymer_precision += rymer_seed_precision

    print(f'Total minimizer seeds: {total_minimizer_seeds}')
    print(f'Total rymer seeds: {total_rymer_seeds}')
    print(f'Minimizer precision: {minimizer_precision / len(reads) if len(reads) > 0 else 0}')
    print(f'Rymer precision: {rymer_precision / len(reads) if len(reads) > 0 else 0}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=1000, help='Number of reads')
    parser.add_argument('--L', type=int, default=100, help='Read length')
    parser.add_argument('--delta', type=float, default=0.01, help='Mutation rate')
    parser.add_argument('--k', type=int, default=15, help='Length of k-mer')
    parser.add_argument('--W', type=int, default=30, help='Length of window')
    args = parser.parse_args()

    assert 0 <= args.delta <= 1, "Mutation rate must be between 0 and 1"
    assert args.N > 0, "Number of reads must be greater than 0"
    assert args.L > 0, "Read length must be greater than 0"
    assert args.k > 0, "Length of k-mer must be greater than 0"
    assert args.W >= args.k, "Length of window must be greater than or equal to length of k-mer"

    main(args.N, args.L, args.delta, args.k, args.W)

