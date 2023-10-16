import os
import re
import csv

def parse_stat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    stats = {}
    for line in lines:
        parts = line.strip().split(":")
        if len(parts) != 2:
            continue
        key, value = parts
        if not value.strip().isdigit():  # Skip lines without a numeric value after the colon
            continue
        stats[key.strip()] = int(value.strip())
    
    return stats

def get_aligner_name_and_damage_type(filename):
    pattern = re.compile(r'.*n([0-9]+)_l([0-9]+)_d([^_]+)_s([\d\.]+)_([\w]+)\.stat')
    match = pattern.match(filename)
    if match:
        numtS_gen, fragment_length, damage_type, subsampling_rate, aligner_name = match.groups()
        return aligner_name, damage_type, fragment_length, subsampling_rate
    else:
        raise ValueError(f'Unexpected file name format: {filename}')

def compute_proportion(directory, output_csv_path):
    files = [f for f in os.listdir(directory) if f.endswith('.stat')]
    
    data = []
    
    for file in files:
        file_path = os.path.join(directory, file)
        stats = parse_stat_file(file_path)
        
        tp = stats.get('True Positives (TP)', 0)
        fp = stats.get('False Positives (FP)', 0)
        tn = stats.get('True Negatives (TN)', 0)
        fn = stats.get('False Negatives (FN)', 0)
        
        tp_mq30 = stats.get('True Positives (TP_MQ30)', 0)
        fp_mq30 = stats.get('False Positives (FP_MQ30)', 0)
        tn_mq30 = stats.get('True Negatives (TN_MQ30)', 0)
        fn_mq30 = stats.get('False Negatives (FN_MQ30)', 0)
        
        aligner_name, damage_type, fragment_length, subsampling_rate = get_aligner_name_and_damage_type(file)
        
        data.append([damage_type, aligner_name, fragment_length, subsampling_rate, tp, fp, tn, fn, tp_mq30, fp_mq30, tn_mq30, fn_mq30])
    
    with open(output_csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        csvwriter.writerow(['Damage_Type', 'Aligner_Name', 'Fragment_Length', 'Subsampling_Rate', 'TP', 'FP', 'TN', 'FN', 'TP_MQ30', 'FP_MQ30', 'TN_MQ30', 'FN_MQ30'])
        
        for row in data:
            csvwriter.writerow(row)

# Define the directory path and output CSV path
directory_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/alignments'
output_csv_path = 'alignment_stats.csv'

# Compute and save the new statistics
compute_proportion(directory_path, output_csv_path)

