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

def get_aligner_name_and_damage_type(filename):
    # Define a regex pattern to capture the damage type and aligner name
    pattern = re.compile(r'.*_l[0-9]+_([^_]+)_s[^_]+_([a-zA-Z0-9]+)\.stat')
    match = pattern.match(filename)
    if match:
        damage_type, aligner_name = match.groups()
        return aligner_name, damage_type
    else:
        raise ValueError(f'Unexpected file name format: {filename}')

def compute_proportion(directory):
    files = [f for f in os.listdir(directory) if f.endswith('.stat')]
    
    # Dictionaries to hold the sum of proportions and count of files by damage type and aligner
    damage_type_aligner_proportion_correct_sum = defaultdict(lambda: defaultdict(float))
    damage_type_aligner_proportion_incorrect_sum = defaultdict(lambda: defaultdict(float))
    damage_type_aligner_file_count = defaultdict(lambda: defaultdict(int))
    
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
        aligner_name, damage_type = get_aligner_name_and_damage_type(file)
        
        # Compute the proportions for this file
        proportion_correct = correct_map / 2000
        proportion_incorrect = 0 if mapped == 0 else (mapped - correct_map) / mapped
        
        # Aggregate the proportions and count of files by damage type and aligner
        damage_type_aligner_proportion_correct_sum[damage_type][aligner_name] += proportion_correct
        damage_type_aligner_proportion_incorrect_sum[damage_type][aligner_name] += proportion_incorrect
        damage_type_aligner_file_count[damage_type][aligner_name] += 1
    
    # Now compute and print the average proportions for each damage type and aligner
    for damage_type in sorted(damage_type_aligner_proportion_correct_sum.keys()):
        print(f'Damage Type: {damage_type}')
        aligner_dict = damage_type_aligner_proportion_correct_sum[damage_type]
        for aligner_name, proportion_correct_sum in aligner_dict.items():
            average_proportion_correct = proportion_correct_sum / damage_type_aligner_file_count[damage_type][aligner_name]
            average_proportion_incorrect = damage_type_aligner_proportion_incorrect_sum[damage_type][aligner_name] / damage_type_aligner_file_count[damage_type][aligner_name]
            print(f'\tAligner: {aligner_name}, Average Proportion of Correctly Mapped Reads: {average_proportion_correct}, Average Proportion of Incorrectly Mapped Reads: {average_proportion_incorrect}')

# Define the directory path
directory_path = '/home/projects/MAAG/Magpie/Magpie/bwa_experiment/human_mito/alignments'

# Compute and print the proportions
compute_proportion(directory_path)

