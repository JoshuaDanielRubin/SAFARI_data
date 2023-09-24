import os
from collections import defaultdict

def parse_stat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    stats = {}
    for line in lines:
        key, value = line.strip().split(":")
        stats[key.strip()] = value.strip()
    
    return stats

def get_aligner_name(filename):
    # Split the filename by underscore and take the last part
    aligner_name_with_extension = filename.split('_')[-1]
    # Now remove the extension to get the aligner name
    aligner_name = aligner_name_with_extension.split('.')[0]
    return aligner_name

def compute_proportion(directory):
    files = [f for f in os.listdir(directory) if f.endswith('.stat')]
    
    # Dictionaries to hold the sum of proportions and count of files by aligner
    aligner_proportion_correct_sum = defaultdict(float)
    aligner_proportion_incorrect_sum = defaultdict(float)
    aligner_file_count = defaultdict(int)
    
    for file in files:
        file_path = os.path.join(directory, file)
        stats = parse_stat_file(file_path)
        correct_map = int(stats['correctmap'])
        mapped = int(stats['mapped'])
        aligner_name = get_aligner_name(file)
        # Compute the proportions for this file
        proportion_correct = correct_map / 10000
        proportion_incorrect = 0 if mapped == 0 else (mapped - correct_map) / mapped
        # Aggregate the proportions and count of files by aligner
        aligner_proportion_correct_sum[aligner_name] += proportion_correct
        aligner_proportion_incorrect_sum[aligner_name] += proportion_incorrect
        aligner_file_count[aligner_name] += 1
    
    # Now compute and print the average proportions for each aligner
    for aligner_name in aligner_proportion_correct_sum.keys():
        average_proportion_correct = aligner_proportion_correct_sum[aligner_name] / aligner_file_count[aligner_name]
        average_proportion_incorrect = aligner_proportion_incorrect_sum[aligner_name] / aligner_file_count[aligner_name]
        print(f'Aligner: {aligner_name}, Average Proportion of Correctly Mapped Reads: {average_proportion_correct}, Average Proportion of Incorrectly Mapped Reads: {average_proportion_incorrect}')

# Define the directory path
directory_path = '/home/projects/MAAG/Magpie/Magpie/bwa_experiment/gabriel_snakefile/alignments'

# Compute and print the proportions
compute_proportion(directory_path)

