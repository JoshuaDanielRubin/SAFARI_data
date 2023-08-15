
import random
import numpy as np
import matplotlib.pyplot as plt

# Function Definitions

def generate_random_dna(length=50000):
    return ''.join(random.choice('ATCG') for _ in range(length))

def create_minimizer_index(dna_sequence, k=20, w=11):
    minimizer_index = {}
    for i in range(len(dna_sequence) - k - w + 2):
        window = dna_sequence[i:i+k+w-1]
        minimizers = [window[j:j+k] for j in range(w)]
        min_kmer = min(minimizers)
        if min_kmer not in minimizer_index:
            minimizer_index[min_kmer] = []
        minimizer_index[min_kmer].append(i + minimizers.index(min_kmer))
    return minimizer_index

def simulate_reads_with_positions(dna_sequence, read_length=100, num_reads=1000):
    reads = []
    positions = []
    for _ in range(num_reads):
        start = random.randint(0, len(dna_sequence) - read_length)
        positions.append(start)
        reads.append(dna_sequence[start:start+read_length])
    return reads, positions

def introduce_deamination(read, rate=0.4):
    mutated_read = []
    for base in read:
        rand_val = random.random()
        if base == "C" and rand_val < rate:
            mutated_read.append("T")
        elif base == "G" and rand_val < rate:
            mutated_read.append("A")
        else:
            mutated_read.append(base)
    return ''.join(mutated_read)

def convert_to_rymer(sequence):
    rymer_sequence = []
    for base in sequence:
        if base in ["C", "T"]:
            rymer_sequence.append("C")
        elif base in ["G", "A"]:
            rymer_sequence.append("A")
        else:
            rymer_sequence.append(base)
    return ''.join(rymer_sequence)

def compute_mismatches(read_kmer, ref_kmer):
    return sum(1 for a, b in zip(read_kmer, ref_kmer) if a != b)

def compute_average_mismatches_for_k(k, w):
    minimizer_index = create_minimizer_index(reference, k, w)
    rymer_reference = convert_to_rymer(reference)
    rymer_index = create_minimizer_index(rymer_reference, k, w)
    
    average_mismatches_list_delta = []
    for delta in delta_values:
        deaminated_reference = introduce_deamination(reference, delta)
        reads, positions = simulate_reads_with_positions(deaminated_reference)
        mismatch_counts_delta = []
        for i, read in enumerate(reads):
            for j in range(len(read) - k + 1):
                kmer = read[j:j+k]
                rymer_kmer = convert_to_rymer(kmer)
                if rymer_kmer in rymer_index and kmer not in minimizer_index:
                    ref_positions = rymer_index[rymer_kmer]
                    for pos in ref_positions:
                        expected_position = positions[i] + j
                        if expected_position == pos:
                            ref_kmer = reference[pos:pos+k]
                            mismatches = compute_mismatches(kmer, ref_kmer)
                            mismatch_counts_delta.append(mismatches)
                            break
                            
        avg_mismatches_delta = np.mean(mismatch_counts_delta) if mismatch_counts_delta else 0
        average_mismatches_list_delta.append(avg_mismatches_delta)
        
    return average_mismatches_list_delta

# Main Script
reference = generate_random_dna()
delta_values = np.arange(0, 1.2, 0.1)
k_values = list(range(5, 21, 2))

plt.figure(figsize=(10, 6))
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange']
labels = []

for idx, k in enumerate(k_values):
    print("k = " + str(k))
    w = k - 2
    average_mismatches_list_delta = compute_average_mismatches_for_k(k, w)
    label = f'k = {k}, w = {w}'
    plt.plot(delta_values, average_mismatches_list_delta, marker='o', linestyle='-', color=colors[idx % len(colors)], label=label)

plt.xlabel('Delta (Deamination Rate)')
plt.ylabel('Average Number of Mismatches')
plt.title('Average Mismatches vs. Delta for Different k values')
plt.legend()
plt.ylim(0, 6)
plt.grid(True)
plt.savefig("avg_mismatches.png")
plt.close()

