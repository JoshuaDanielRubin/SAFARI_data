import os
from collections import defaultdict
import re

def parse_stat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    stats = {}
    for line in lines:
        key, value = line.strip().split(":")
        stats[key.strip()] = value.strip()
    
    return stats

def get_aligner_name_and_sampling_rate(filename):
    # Define a regex pattern to capture the sampling rate and aligner name
    pattern = re.compile(r'.*_s([^_]+)_([^\.]+)\.stat')
    match = pattern.match(filename)
    if match:
        sampling_rate, aligner_name = match.groups()
        return aligner_name, sampling_rate
    else:
        raise ValueError(f'Unexpected file name format: {filename}')

def compute_proportion(directory):
    files = [f for f in os.listdir(directory) if f.endswith('.stat')]
    
    # Dictionaries to hold the sum of proportions and count of files by sampling rate and aligner
    sampling_aligner_proportion_correct_sum = defaultdict(lambda: defaultdict(float))
    sampling_aligner_proportion_incorrect_sum = defaultdict(lambda: defaultdict(float))
    sampling_aligner_file_count = defaultdict(lambda: defaultdict(int))
    
    for file in files:
        file_path = os.path.join(directory, file)
        stats = parse_stat_file(file_path)
        total = stats.get('Total')
        if total is not None:
            total = int(total)
        else:
            continue
        correct_map = int(stats['correctmap'])
        mapped = int(stats['mapped'])
        aligner_name, sampling_rate = get_aligner_name_and_sampling_rate(file)
        # Compute the proportions for this file
        proportion_correct = correct_map / 1000
        proportion_incorrect = 0 if mapped == 0 else (mapped - correct_map) / mapped
        # Aggregate the proportions and count of files by sampling rate and aligner
        sampling_aligner_proportion_correct_sum[sampling_rate][aligner_name] += proportion_correct
        sampling_aligner_proportion_incorrect_sum[sampling_rate][aligner_name] += proportion_incorrect
        sampling_aligner_file_count[sampling_rate][aligner_name] += 1
    
    # Now compute and print the average proportions for each sampling rate and aligner
    for sampling_rate in sorted(sampling_aligner_proportion_correct_sum.keys(), key=lambda x: float(x)):
        print(f'Sampling Rate: {sampling_rate}')
        aligner_dict = sampling_aligner_proportion_correct_sum[sampling_rate]
        for aligner_name, proportion_correct_sum in aligner_dict.items():
            average_proportion_correct = proportion_correct_sum / sampling_aligner_file_count[sampling_rate][aligner_name]
            average_proportion_incorrect = sampling_aligner_proportion_incorrect_sum[sampling_rate][aligner_name] / sampling_aligner_file_count[sampling_rate][aligner_name]
            print(f'\tAligner: {aligner_name}, Average Proportion of Correctly Mapped Reads: {average_proportion_correct}, Average Proportion of Incorrectly Mapped Reads: {average_proportion_incorrect}')

# Define the directory path
directory_path = '/home/projects/MAAG/Magpie/Magpie/bwa_experiment/gabriel_snakefile/alignments'

# Compute and print the proportions
compute_proportion(directory_path)

