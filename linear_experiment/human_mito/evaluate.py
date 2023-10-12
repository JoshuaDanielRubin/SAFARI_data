import os
import re
import csv

def parse_stat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    stats = {}
    for line in lines:
        if line.startswith('-----------------------------------------'):
            continue
        parts = line.strip().split(":")
        if len(parts) != 2:
            continue
        key, value = parts
        stats[key.strip()] = int(value.strip())
    
    return stats

def get_aligner_name_and_damage_type(filename):
    # Updated pattern for the new naming format
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
        
        total_reads = stats.get('Total reads', 0)
        unmapped_reads = stats.get('Unmapped reads', 0)  # Extracting the unmapped reads
        mapped_to_mt = stats.get('Mapped to MT', 0)
        mapped_to_mt_correct_location = stats.get('Mapped to MT (Correct Location)', 0)
        mapped_to_mt_correct_location_mq_30 = stats.get('Mapped to MT (Correct Location, MQ>30)', 0)
        mapped_not_to_mt = stats.get('NOT Mapped to MT', 0)
        mapped_not_to_mt_mq_30 = stats.get('NOT Mapped to MT (MQ>30)', 0)
        
        aligner_name, damage_type, fragment_length, subsampling_rate = get_aligner_name_and_damage_type(file)
        
        data.append([damage_type, aligner_name, fragment_length, subsampling_rate, total_reads, unmapped_reads, mapped_to_mt, 
                     mapped_to_mt_correct_location, mapped_to_mt_correct_location_mq_30, mapped_not_to_mt, mapped_not_to_mt_mq_30])
    
    with open(output_csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        csvwriter.writerow(['Damage_Type', 'Aligner_Name', 'Fragment_Length', 'Subsampling_Rate', 'Total_Reads', 'Unmapped_Reads', 'Mapped_to_MT', 
                            'Mapped_to_MT_Correct_Location', 'Mapped_to_MT_Correct_Location_MQ>30', 'Mapped_NOT_to_MT', 'Mapped_NOT_to_MT_MQ>30'])
        
        for row in data:
            csvwriter.writerow(row)

# Define the directory path and output CSV path
directory_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/alignments'
output_csv_path = 'alignment_stats.csv'

# Compute and save the new statistics
compute_proportion(directory_path, output_csv_path)

