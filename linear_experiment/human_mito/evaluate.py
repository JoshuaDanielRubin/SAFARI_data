import os
from collections import defaultdict
import re
import csv

def parse_stat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    stats = {}
    for line in lines:
        key, value = line.strip().split(":")
        stats[key.strip()] = value.strip()
    
    return stats

def get_aligner_name_and_damage_type(filename):
    # Modified regex pattern to capture aligner names with underscores
    pattern = re.compile(r'.*_l[0-9]+_([^_]+)_s[^_]+_([\w]+)\.stat')
    match = pattern.match(filename)
    if match:
        damage_type, aligner_name = match.groups()
        return aligner_name, damage_type
    else:
        raise ValueError(f'Unexpected file name format: {filename}')

def compute_proportion(directory, output_csv_path):
    files = [f for f in os.listdir(directory) if f.endswith('.stat')]
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
        proportion_correct = correct_map / 1000
        proportion_incorrect = 0 if mapped == 0 else (mapped - correct_map) / mapped
        damage_type_aligner_proportion_correct_sum[damage_type][aligner_name] += proportion_correct
        damage_type_aligner_proportion_incorrect_sum[damage_type][aligner_name] += proportion_incorrect
        damage_type_aligner_file_count[damage_type][aligner_name] += 1
    
    # Save the results to a CSV file
    with open(output_csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Damage_Type', 'Aligner_Name', 'Avg_Proportion_Correct', 'Avg_Proportion_Incorrect'])
        
        for damage_type in sorted(damage_type_aligner_proportion_correct_sum.keys()):
            aligner_dict = damage_type_aligner_proportion_correct_sum[damage_type]
            for aligner_name, proportion_correct_sum in aligner_dict.items():
                average_proportion_correct = proportion_correct_sum / damage_type_aligner_file_count[damage_type][aligner_name]
                average_proportion_incorrect = damage_type_aligner_proportion_incorrect_sum[damage_type][aligner_name] / damage_type_aligner_file_count[damage_type][aligner_name]
                csvwriter.writerow([damage_type, aligner_name, average_proportion_correct, average_proportion_incorrect])

# Define the directory path and output CSV path
directory_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/alignments'
output_csv_path = 'average_proportions.csv'

# Compute and save the proportions
compute_proportion(directory_path, output_csv_path)

