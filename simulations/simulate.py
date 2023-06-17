import random
from collections import defaultdict
import argparse
import csv
import os

def print_stats(params, stats):
    labels = ['N', 'L', 'delta', 'k', 'W', 'Total minimizer seeds', 'Total rymer seeds',
              'Minimizer precision', 'Rymer precision', 'Rymer spuriousness',
              'Rymer recovery rate', 'Injectivity']
    for label, value in zip(labels, params + stats):
        print(f'{label}: {value}')

def generate_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def reverse_complement(kmer):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(kmer))

def create_minimizer_index(S, W, k):
    minimizer_index = defaultdict(list)
    rymer_index = defaultdict(list)
    rymer_minimizer_map = {}

    for i in range(len(S) - W + 1):
        window = S[i:i+W]
        for j in range(len(window) - k + 1):
            kmer = window[j:j+k]
            minimizer = min(kmer, reverse_complement(kmer))
            rymer = ''.join(['R' if base in 'GA' else 'Y' for base in kmer])
            minimizer_index[minimizer].append(i + j)
            rymer_index[rymer].append(i + j)
            if rymer not in rymer_minimizer_map:
                rymer_minimizer_map[rymer] = set()
            rymer_minimizer_map[rymer].add(minimizer)

    return minimizer_index, rymer_index, rymer_minimizer_map

def generate_reads(S, N, L, delta):
    reads = []
    for _ in range(N):
        start = random.randint(0, len(S) - L)
        read = list(S[start:start+L])
        deaminated_bases = set()
        for i in range(len(read)):
            if (read[i] == 'C' or read[i] == 'G') and random.random() < delta:
                if read[i] == 'C':
                    read[i] = 'T'
                elif read[i] == 'G':
                    read[i] = 'A'
                deaminated_bases.add(i)
        reads.append((''.join(read), start, deaminated_bases))
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

def calculate_injectivity(rymer_minimizer_map):
    unique_minimizers = set()
    total_rymers = len(rymer_minimizer_map)
    for minimizers in rymer_minimizer_map.values():
        unique_minimizers.update(minimizers)
    injectivity = total_rymers / len(unique_minimizers)
    return injectivity

def compute_rymer_spuriousness(minimizer_seeds, rymer_seeds):
    minimizer_kmers = {pos for pos, _ in minimizer_seeds}
    rymer_kmers = {pos for pos, correct in rymer_seeds if correct}

    # Kmers where the minimizer does not match but the rymer does, and is correct
    spurious_rymer_kmers = rymer_kmers.difference(minimizer_kmers)

    # Fraction of these kmers over the total number of rymer seeds
    if rymer_kmers:
        rymer_spuriousness = 1 - (len(spurious_rymer_kmers) / len(rymer_kmers))
    else:
        rymer_spuriousness = 0

    return rymer_spuriousness

def compute_rymer_recovery(minimizer_seeds, rymer_seeds, deaminated_bases, k):
    minimizer_kmers = {pos for pos, _ in minimizer_seeds}
    rymer_kmers = {pos for pos, correct in rymer_seeds if correct}

    # Kmers where the minimizer does not match but the rymer does, and is correct
    spurious_rymer_kmers = rymer_kmers.difference(minimizer_kmers)

    # Of these, the ones that can be explained by deamination
    deaminated_rymer_kmers = {pos for pos in spurious_rymer_kmers for i in range(k) if pos+i in deaminated_bases}

    # Fraction of these kmers over the total number of spurious rymer kmers
    if spurious_rymer_kmers:
        rymer_recovery_rate = len(deaminated_rymer_kmers) / len(spurious_rymer_kmers)
    else:
        rymer_recovery_rate = 0

    return rymer_recovery_rate

def write_header(filename):
    header = ['N', 'L', 'delta', 'k', 'W', 'Total minimizer seeds', 'Total rymer seeds',
              'Minimizer precision', 'Rymer precision', 'Rymer spuriousness',
              'Rymer recovery rate', 'Injectivity']

    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        with open(filename, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(header)

def append_row(filename, params, stats):
    row = [*params, *stats]
    with open(filename, 'a', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(row)

def main(N, L, delta, k, W, stdout):
    S = generate_dna(10000)
    minimizer_index, rymer_index, rymer_minimizer_map = create_minimizer_index(S, W, k)
    reads = generate_reads(S, N, L, delta)

    total_minimizer_seeds = 0
    total_rymer_seeds = 0
    minimizer_precision = 0
    rymer_precision = 0
    rymer_spuriousness_sum = 0
    rymer_recovery_sum = 0

    for read, origin, deaminated_bases in reads:
        minimizer_seeds = find_seeds(minimizer_index, read, k, origin)
        rymer_read = ''.join(['R' if base in 'GA' else 'Y' for base in read])
        rymer_seeds = find_seeds(rymer_index, rymer_read, k, origin)

        minimizer_seeds_found, minimizer_seed_precision = compute_stats(minimizer_seeds)
        rymer_seeds_found, rymer_seed_precision = compute_stats(rymer_seeds)

        total_minimizer_seeds += minimizer_seeds_found
        total_rymer_seeds += rymer_seeds_found
        minimizer_precision += minimizer_seed_precision
        rymer_precision += rymer_seed_precision
        rymer_spuriousness = compute_rymer_spuriousness(minimizer_seeds, rymer_seeds)
        rymer_spuriousness_sum += rymer_spuriousness
        rymer_recovery_rate = compute_rymer_recovery(minimizer_seeds, rymer_seeds, deaminated_bases, k)
        rymer_recovery_sum += rymer_recovery_rate

    injectivity = calculate_injectivity(rymer_minimizer_map)

    params = [N, L, delta, k, W]
    stats = [
    total_minimizer_seeds,
    total_rymer_seeds,
    "{:.3f}".format(minimizer_precision / len(reads)) if len(reads) > 0 else 0,
    "{:.3f}".format(rymer_precision / len(reads)) if len(reads) > 0 else 0,
    "{:.3f}".format(rymer_spuriousness_sum / len(reads)) if len(reads) > 0 else 0,
    "{:.3f}".format(rymer_recovery_sum / len(reads)) if len(reads) > 0 else 0,
    "{:.3f}".format(injectivity)
            ]

    if stdout:
        print_stats(params, stats)
    else:
        filename = 'results.tsv'
        write_header(filename)
        append_row(filename, params, stats)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=1000, help='Number of reads')
    parser.add_argument('--L', type=int, default=100, help='Read length')
    parser.add_argument('--delta', type=float, default=0.01, help='Mutation rate')
    parser.add_argument('--k', type=int, default=15, help='Length of k-mer')
    parser.add_argument('--W', type=int, default=30, help='Length of window')
    parser.add_argument('--stdout', action='store_true', help='Print results to stdout instead of writing to a file')
    args = parser.parse_args()

    assert 0 <= args.delta <= 1, "Mutation rate must be between 0 and 1"
    assert args.N > 0, "Number of reads must be greater than 0"
    assert args.L > 0, "Read length must be greater than 0"
    assert args.k > 0, "Length of k-mer must be greater than 0"
    assert args.W >= args.k, "Length of window must be greater than or equal to length of k-mer"

    main(args.N, args.L, args.delta, args.k, args.W, args.stdout)

