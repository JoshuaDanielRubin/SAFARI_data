import os
import glob
import pandas as pd
import numpy as np
import io
import re
import matplotlib.pyplot as plt
from collections import defaultdict

def clean_data(df):
    assert df is not None, "DataFrame should not be None."
    df = df.applymap(lambda x: re.search(r"([\d\.]+)", str(x)).group(1) if re.search(r"([\d\.]+)", str(x)) else np.nan)
    return df.astype(float)

def extract_damage_type(file_name):
    parts = file_name.split('_')
    for part in parts:
        if part.startswith('d') or part.startswith('dd'):
            damage_type = part.lstrip('d').lstrip('d')
            assert damage_type in ['none', 'mid', 'high', 'single'], f'Unexpected damage type: {damage_type}'
            return damage_type
    raise ValueError(f'Unexpected file name structure, no damage type found: {file_name}')

def extract_aligner(file_name):
    return file_name.split('_')[-1].split('.')[0]

def load_damage_data(file_name):
    file_path = os.path.join(damage_data_path, file_name)
    assert os.path.exists(file_path), f"File path does not exist: {file_path}"
    if os.path.getsize(file_path) == 0:
        print(f'File {file_name} is empty.')
        return None
    try:
        data = pd.read_csv(file_path, delimiter='\t', index_col=0)
        assert not data.empty, f"No data in file {file_name}."
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None
    return data

def load_prof_data(file_name):
    file_path = os.path.join(prof_data_path, file_name)
    assert os.path.exists(file_path), f"File path does not exist: {file_path}"
    if os.path.getsize(file_path) == 0:
        return None, None
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            second_header_index = len(lines) // 2
            assert second_header_index > 0, "Invalid file format."
            table1_lines = lines[:second_header_index]
            table2_lines = lines[second_header_index:]

        table1 = pd.read_csv(io.StringIO(''.join(table1_lines)), delimiter='\t', header=0)
        table2 = pd.read_csv(io.StringIO(''.join(table2_lines)), delimiter='\t', header=0)
        assert not (table1.empty or table2.empty), f"No data in file {file_name}."
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None, None
    return table1, table2

def normalize_data(df):
    row_sums = df.sum(axis=1)
    return df.div(row_sums, axis=0)

def compute_kl_divergence(true_data, estimated_data, aligner, damage_type):
    assert true_data is not None, 'True data is missing.'
    assert estimated_data is not None, 'Estimated data is missing.'
    
    if np.all(true_data.isna()) or np.all(estimated_data.isna()):
        return None, 0
    
    true_data = clean_data(true_data)
    estimated_data = clean_data(estimated_data)

    # Normalize the data
    true_data = normalize_data(true_data)
    estimated_data = normalize_data(estimated_data)
    
    assert true_data.shape == estimated_data.shape, "Shape mismatch between true and estimated data."

    nan_indices = true_data.isnull().any(axis=1) | estimated_data.isnull().any(axis=1)
    true_data = true_data[~nan_indices]
    estimated_data = estimated_data[~nan_indices]

    if true_data.empty or estimated_data.empty:
        return None, 0

    epsilon = 1e-9  
    flat_true_data = true_data.values.flatten() + epsilon  
    flat_estimated_data = estimated_data.values.flatten() + epsilon  

    kl_divergence = np.mean(flat_true_data * np.log(flat_true_data / flat_estimated_data))
    
    if kl_divergence < 0:
        print(f"Negative KL divergence detected: {kl_divergence}")
        print(f"True data: {flat_true_data}")
        print(f"Estimated data: {flat_estimated_data}")
        return None, 0  # or you can raise an exception if you prefer

    return kl_divergence, len(true_data)

def filter_kl_data_for_giraffe_and_safari(kl_data):
    filtered_kl_data = defaultdict(lambda: defaultdict(float))
    for aligner in ['giraffe', 'safari']:
        if aligner in kl_data:
            filtered_kl_data[aligner] = kl_data[aligner]
    return filtered_kl_data

def plot_kl(kl_data, kl_sample_count_data, plot_title, save_file_name):
    colorblind_colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', '#CA9161', '#FBAFE4']
    aligners = list(kl_data.keys())
    damage_types = list(kl_data[aligners[0]].keys())
    
    x = np.arange(len(aligners))
    width = 0.2

    fig, ax = plt.subplots()
    for i, damage_type in enumerate(damage_types):
        kl_values = [kl_data[aligner][damage_type] for aligner in aligners]
        ax.bar(x + i*width, kl_values, width, label=damage_type, color=colorblind_colors[i])
    
    ax.set_xlabel('Aligner')
    ax.set_ylabel('Average KL Divergence')
    ax.set_title(plot_title)
    ax.set_xticks(x + width*(len(damage_types)-1)/2)
    ax.set_xticklabels(aligners)
    ax.legend()

    fig.tight_layout()
    plt.savefig(save_file_name)
    
    for aligner in aligners:
        for damage_type in damage_types:
            sample_count = kl_sample_count_data[aligner][damage_type]
            print(f'Average kl for {aligner} with damage_type {damage_type}: {kl_data[aligner][damage_type]} (based on {sample_count} samples)')

# The `check_data` function remains unchanged
def check_data(damage_data_dict, prof_data_dict):
    kl_sum_data = defaultdict(lambda: defaultdict(float))
    kl_count_data = defaultdict(lambda: defaultdict(int))
    kl_sample_count_data = defaultdict(lambda: defaultdict(int))
    for file_name, (table1, table2) in prof_data_dict.items():
        damage_type = extract_damage_type(file_name)
        aligner = extract_aligner(file_name)

        if table1 is not None:
            true_data_key = f'{damage_type}{len(table1)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            kl_table1, sample_count_table1 = compute_kl_divergence(true_data, table1, aligner, damage_type)

            if kl_table1 is not None:
                kl_sum_data[aligner][damage_type] += kl_table1
                kl_count_data[aligner][damage_type] += 1
                kl_sample_count_data[aligner][damage_type] += sample_count_table1

        if table2 is not None:
            true_data_key = f'{damage_type}{len(table2)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            kl_table2, sample_count_table2 = compute_kl_divergence(true_data, table2, aligner, damage_type)

            if kl_table2 is not None:
                kl_sum_data[aligner][damage_type] += kl_table2
                kl_count_data[aligner][damage_type] += 1
                kl_sample_count_data[aligner][damage_type] += sample_count_table2

    kl_avg_data = defaultdict(lambda: defaultdict(float))
    for aligner, damage_data in kl_sum_data.items():
        for damage_type, kl_sum in damage_data.items():
            kl_avg_data[aligner][damage_type] = kl_sum / kl_count_data[aligner][damage_type]

    #plot_kl(kl_avg_data, kl_sample_count_data, 'Average KL Divergence by Aligner and Damage Type', 'kl_plot.png')
    filtered_kl_data = filter_kl_data_for_giraffe_and_safari(kl_avg_data)
    plot_kl(filtered_kl_data, kl_sample_count_data, 'Average KL Divergence between Giraffe and SAFARI', 'kl_plot_giraffe_safari.png')

if __name__ == "__main__":
    damage_data_path = '.'  # Your path to damage data files
    prof_data_path = 'new_alignments/profs'  # Your path to prof data files

    assert os.path.exists(damage_data_path), "Damage data path does not exist."
    assert os.path.exists(prof_data_path), "Prof data path does not exist."

    damage_data_files = glob.glob(os.path.join(damage_data_path, '*.dat'))
    prof_data_files = glob.glob(os.path.join(prof_data_path, '*.prof'))

    assert damage_data_files, "No damage data files found."
    assert prof_data_files, "No prof data files found."

    damage_data_dict = {os.path.basename(file_name): load_damage_data(os.path.basename(file_name)) for file_name in damage_data_files}
    prof_data_dict = {os.path.basename(file_name): load_prof_data(os.path.basename(file_name)) for file_name in prof_data_files}

    check_data(damage_data_dict, prof_data_dict)
