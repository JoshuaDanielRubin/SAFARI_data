import argparse
import csv
import os
import random

from minimizer_sketch import MinimizerSketch
from read_generator import ReadGenerator
from stats_calculator import StatsCalculator

def generate_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def print_stats(params, stats):
    labels = ['N', 'L', 'delta', 'k', 'W', 'Total minimizer seeds', 'Total rymer seeds',
              'Minimizer precision', 'Rymer precision', 'Rymer spuriousness',
              'Rymer recovery rate', 'Injectivity']
    for label, value in zip(labels, params + stats):
        print(f'{label}: {value}')


def calculate_injectivity(rymer_minimizer_map):
    unique_minimizers = set()
    total_rymers = len(rymer_minimizer_map)
    for minimizers in rymer_minimizer_map.values():
        unique_minimizers.update(minimizers)
    injectivity = total_rymers / len(unique_minimizers)
    return injectivity


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
    sketch = MinimizerSketch(S, W, k)
    read_generator = ReadGenerator(S, delta)
    reads = read_generator.generate_reads(N, L)

    total_minimizer_seeds = 0
    total_rymer_seeds = 0
    minimizer_precision = 0
    rymer_precision = 0
    rymer_spuriousness_sum = 0
    rymer_recovery_sum = 0

    for read, origin, deaminated_bases in reads:
        minimizer_seeds = sketch.find_seeds(sketch.minimizer_index, read, origin)
        rymer_read = ''.join(['R' if base in 'GA' else 'Y' for base in read])
        rymer_seeds = sketch.find_seeds(sketch.rymer_index, rymer_read, origin)

        stats_calculator = StatsCalculator(minimizer_seeds, rymer_seeds, deaminated_bases, k)
        total_minimizer, total_rymer, minimizer_prec, rymer_prec, rymer_spuriousness, rymer_recovery = stats_calculator.compute_stats()

        total_minimizer_seeds += total_minimizer
        total_rymer_seeds += total_rymer
        minimizer_precision += minimizer_prec
        rymer_precision += rymer_prec
        rymer_spuriousness_sum += rymer_spuriousness
        rymer_recovery_sum += rymer_recovery

    injectivity = calculate_injectivity(sketch.rymer_minimizer_map)

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

