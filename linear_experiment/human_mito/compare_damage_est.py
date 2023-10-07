import os
import glob
import pandas as pd
import numpy as np
import io
import re
import matplotlib.pyplot as plt
from collections import defaultdict
from math import sqrt
from statistics import median
from sklearn.metrics.pairwise import cosine_similarity

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

def compute_cosine_similarity(true_data, estimated_data, aligner, damage_type):
    assert true_data is not None, 'True data is missing.'
    assert estimated_data is not None, 'Estimated data is missing.'

    if np.all(true_data.isna()) or np.all(estimated_data.isna()):
        return None, 0

    true_data = clean_data(true_data)
    estimated_data = clean_data(estimated_data)

    assert true_data.shape == estimated_data.shape, "Shape mismatch between true and estimated data."

    nan_indices = true_data.isnull().any(axis=1) | estimated_data.isnull().any(axis=1)
    true_data = true_data[~nan_indices]
    estimated_data = estimated_data[~nan_indices]

    if true_data.empty or estimated_data.empty:
        return None, 0

    # Handle rows with all zeros explicitly
    true_data_zero_rows = np.all(true_data == 0, axis=1)
    estimated_data_zero_rows = np.all(estimated_data == 0, axis=1)

    both_zero_rows = true_data_zero_rows & estimated_data_zero_rows
    one_zero_rows = true_data_zero_rows ^ estimated_data_zero_rows

    # Debug prints
    print(f"Number of both zero rows: {np.sum(both_zero_rows)}")
    print(f"Number of one zero rows: {np.sum(one_zero_rows)}")

    non_zero_rows = ~(true_data_zero_rows | estimated_data_zero_rows)

    # Compute cosine similarity only for non-zero rows
    if np.any(non_zero_rows):
        cosine_sim_non_zero = cosine_similarity(true_data[non_zero_rows], estimated_data[non_zero_rows])
    else:
        cosine_sim_non_zero = np.array([])

    # Combine results
    cosine_sim = np.concatenate([
        np.ones(np.sum(both_zero_rows)),
        np.zeros(np.sum(one_zero_rows)),
        cosine_sim_non_zero.flatten()
    ])
    
    mean_cosine_sim = np.mean(cosine_sim)

    return mean_cosine_sim, len(true_data)

def filter_cos_sim_data_for_giraffe_and_safari(cos_sim_data):
    filtered_cos_sim_data = defaultdict(lambda: defaultdict(float))
    for aligner in ['giraffe', 'safari']:
        if aligner in cos_sim_data:
            filtered_cos_sim_data[aligner] = cos_sim_data[aligner]
    return filtered_cos_sim_data

def plot_rmse(cos_sim_data, cos_sim_sample_count_data, plot_title, save_file_name):
    colorblind_colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', '#CA9161', '#FBAFE4']
    aligners = list(cos_sim_data.keys())
    damage_types = list(cos_sim_data[aligners[0]].keys())
    
    x = np.arange(len(aligners))
    width = 0.2

    fig, ax = plt.subplots()
    for i, damage_type in enumerate(damage_types):
        cos_sim_values = [cos_sim_data[aligner][damage_type] for aligner in aligners]
        ax.bar(x + i*width, cos_sim_values, width, label=damage_type, color=colorblind_colors[i])
    
    ax.set_xlabel('Aligner')
    ax.set_ylabel('Average RMSE')
    ax.set_title(plot_title)
    ax.set_xticks(x + width*(len(damage_types)-1)/2)
    ax.set_xticklabels(aligners)
    ax.legend()

    fig.tight_layout()
    plt.savefig(save_file_name)
    
    for aligner in aligners:
        for damage_type in damage_types:
            sample_count = cos_sim_sample_count_data[aligner][damage_type]
            print(f'Average rmse for {aligner} with damage_type {damage_type}: {cos_sim_data[aligner][damage_type]} (based on {sample_count} samples)')

# The `check_data` function remains unchanged
def check_data(damage_data_dict, prof_data_dict):
    cos_sim_values_data = defaultdict(lambda: defaultdict(list))
    cos_sim_sample_count_data = defaultdict(lambda: defaultdict(int))
    
    for file_name, (table1, table2) in prof_data_dict.items():
        damage_type = extract_damage_type(file_name)
        aligner = extract_aligner(file_name)

        if table1 is not None:
            true_data_key = f'{damage_type}{len(table1)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            cos_sim_table1, sample_count_table1 = compute_cosine_similarity(true_data, table1, aligner, damage_type)

            if cos_sim_table1 is not None:
                cos_sim_values_data[aligner][damage_type].append(cos_sim_table1)
                cos_sim_sample_count_data[aligner][damage_type] += sample_count_table1

        if table2 is not None:
            true_data_key = f'{damage_type}{len(table2)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            cos_sim_table2, sample_count_table2 = compute_cosine_similarity(true_data, table2, aligner, damage_type)

            if cos_sim_table2 is not None:
                cos_sim_values_data[aligner][damage_type].append(cos_sim_table2)
                cos_sim_sample_count_data[aligner][damage_type] += sample_count_table2

    cos_sim_median_data = defaultdict(lambda: defaultdict(float))
    for aligner, damage_data in cos_sim_values_data.items():
        for damage_type, cos_sim_values in damage_data.items():
            cos_sim_median_data[aligner][damage_type] = median(cos_sim_values)

    filtered_cos_sim_data = filter_cos_sim_data_for_giraffe_and_safari(cos_sim_median_data)
    plot_rmse(filtered_cos_sim_data, cos_sim_sample_count_data, 'Median RMSE between Giraffe and SAFARI', 'cos_sim_plot_giraffe_safari.png')


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
