import argparse
import csv
import os
import random

from minimizer_sketch import MinimizerSketch, UniqueMinimizerSketch
from read_generator import ReadGenerator
from stats_calculator import StatsCalculator

def generate_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def print_stats(params, stats):
    labels = ['N', 'L', 'delta', 'k', 'W', 'Total minimizer seeds', 'Total rymer seeds',
              'Minimizer precision', 'Rymer precision', 'Incorrect Rescue Rate',
              'Correct Rescue Rate', 'Uniqueness', 'Minimizer sketch']
    for label, value in zip(labels, params + stats):
        print(f'{label}: {value}')

def calculate_uniqueness(rymer_minimizer_map):
    unique_minimizers = set()
    total_rymers = len(rymer_minimizer_map)
    for minimizers in rymer_minimizer_map.values():
        unique_minimizers.update(minimizers)
    uniqueness = total_rymers / len(unique_minimizers)
    return uniqueness

def write_header(filename):
    header = ['N', 'L', 'delta', 'k', 'W', 'Total minimizer seeds', 'Total rymer seeds',
              'Minimizer precision', 'Rymer precision', 'Incorrect Rescue Rate',
              'Correct Rescue Rate', 'Uniqueness', 'Minimizer sketch']

    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        with open(filename, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(header)

def append_row(filename, params, stats):
    row = [*params, *stats]
    with open(filename, 'a', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(row)

def get_reference_dna(filename):
    with open(filename, 'r') as file:
        # Skip the fasta header
        next(file)
        return ''.join(line.strip() for line in file)

def main(N, L, delta, k, W, unique_minimizer_sketch, stdout, contamination_rate):

    if contamination_rate == 0:
        S = get_reference_dna("rCRS.fa")
    else:
        reference_length = int((1 - contamination_rate) * 2000)
        S = get_reference_dna("rCRS.fa")[:reference_length] + generate_dna(2000 - reference_length)
    
    if unique_minimizer_sketch:
        sketch = UniqueMinimizerSketch(S, W, k)
        sketch_type = 'Unique'
    else:
        sketch = MinimizerSketch(S, W, k)
        sketch_type = 'Naive'
    
    read_generator = ReadGenerator(S, delta)
    reads = read_generator.generate_reads(N, L)
    deaminated_bases = read_generator.deaminated_bases  # Collect deaminated_bases

    total_minimizer_seeds = len(sketch.minimizer_index)
    total_rymer_seeds = len(sketch.rymer_index)
    minimizer_precision = 0
    rymer_precision = 0
    rymer_incorrect_rescue_sum = 0
    rymer_correct_rescue_sum = 0

    for read, origin, _ in reads:
        minimizer_seeds = sketch.find_seeds(sketch.minimizer_index, read, origin)
        rymer_read = ''.join(['R' if base in 'GA' else 'Y' for base in read])
        rymer_seeds = sketch.find_seeds(sketch.rymer_index, rymer_read, origin)

        stats_calculator = StatsCalculator(minimizer_seeds, rymer_seeds, deaminated_bases, k)  # Pass deaminated_bases
        _, _, minimizer_prec, rymer_prec, rymer_incorrect_rescue, rymer_correct_rescue = stats_calculator.compute_stats()

        minimizer_precision += minimizer_prec
        rymer_precision += rymer_prec
        rymer_incorrect_rescue_sum += rymer_incorrect_rescue
        rymer_correct_rescue_sum += rymer_correct_rescue

    uniqueness = calculate_uniqueness(sketch.rymer_minimizer_map)

    params = [N, L, delta, k, W]
    stats = [
        total_minimizer_seeds,
        total_rymer_seeds,
        "{:.3f}".format(minimizer_precision / len(reads)) if len(reads) > 0 else 0,
        "{:.3f}".format(rymer_precision / len(reads)) if len(reads) > 0 else 0,
        "{:.3f}".format(rymer_incorrect_rescue_sum / len(reads)) if len(reads) > 0 else 0,
        "{:.3f}".format(rymer_correct_rescue_sum / len(reads)) if len(reads) > 0 else 0,
        "{:.3f}".format(uniqueness),
        sketch_type
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
    parser.add_argument('--delta', type=float, default=0.3, help='Mutation rate')
    parser.add_argument('--k', type=int, default=15, help='Length of k-mer')
    parser.add_argument('--W', type=int, default=30, help='Length of window')
    parser.add_argument('--unique', action='store_true', help='Use unique minimizer sketch')
    parser.add_argument('--stdout', action='store_true', help='Print results to stdout instead of writing to a file')
    parser.add_argument('--contamination_rate', type=float, default=0.0, help='Proportion of DNA that should be random')
    args = parser.parse_args()

    assert 0 <= args.delta <= 1, "Mutation rate must be between 0 and 1"
    assert args.N > 0, "Number of reads must be greater than 0"
    assert args.L > 0, "Read length must be greater than 0"
    assert args.k > 0, "Length of k-mer must be greater than 0"
    assert args.W >= args.k, "Length of window must be greater than or equal to length of k-mer"

    main(args.N, args.L, args.delta, args.k, args.W, args.unique, args.stdout, args.contamination_rate)
