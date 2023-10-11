import os
import glob
import pandas as pd
import numpy as np
import io
import re
import matplotlib.pyplot as plt
from collections import defaultdict
from math import sqrt

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

def mean_squared_error(true_values, predicted_values):
    assert len(true_values) == len(predicted_values), "Length mismatch between true and predicted values."
    return np.mean((true_values - predicted_values)**2)

def compute_rmse_divergence(true_data, estimated_data, aligner, damage_type):
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

    flat_true_data = true_data.values.flatten()
    flat_estimated_data = estimated_data.values.flatten()

    rmse = sqrt(mean_squared_error(flat_true_data, flat_estimated_data))

    return rmse, len(true_data)

def filter_rmse_data_for_giraffe_and_safari(rmse_data):
    filtered_rmse_data = defaultdict(lambda: defaultdict(float))
    for aligner in ['giraffe', 'safari']:
        if aligner in rmse_data:
            filtered_rmse_data[aligner] = rmse_data[aligner]
    return filtered_rmse_data

def plot_rmse(rmse_data, rmse_sample_count_data, plot_title, save_file_name):
    grey_color = '#808080'
    aligners = ['giraffe', 'safari'] + [a for a in rmse_data.keys() if a not in ['giraffe', 'safari']]
    display_labels = {'giraffe': 'Giraffe', 'safari': 'SAFARI'}
    damage_types = list(rmse_data[aligners[0]].keys())
    
    x = np.arange(len(aligners))
    width = 0.2

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Median RMSE between Estimated and Ground Truth Damage Rate Matrices", fontsize=14)

    for i, damage_type in enumerate(damage_types):
        ax = axs[i // 2, i % 2]
        rmse_values = [rmse_data[aligner][damage_type] for aligner in aligners]
        
        for j, aligner in enumerate(aligners):
            if aligner == 'giraffe':
                color = 'orange'
            elif aligner == 'safari':
                color = 'green'
            else:
                color = grey_color  # Grey color for other aligners
            label = display_labels.get(aligner, aligner)
            ax.bar(x[j], rmse_values[j], width, label=label, color=color)
        
        ax.set_xlabel('Aligner')
        ax.set_ylabel('Median RMSE')
        ax.set_title(f"Damage Matrix: {damage_type.capitalize()}")
        ax.set_xticks(x)
        ax.set_xticklabels([display_labels.get(aligner, aligner) for aligner in aligners], fontweight='bold' if aligner in ['giraffe', 'safari'] else 'normal')

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(save_file_name)
    
    for aligner in aligners:
        for damage_type in damage_types:
            sample_count = rmse_sample_count_data[aligner][damage_type]
            print(f'Median RMSE for {aligner} with damage_type {damage_type}: {rmse_data[aligner][damage_type]} (based on {sample_count} samples)')


def check_data(damage_data_dict, prof_data_dict):
    rmse_values_data = defaultdict(lambda: defaultdict(list))
    rmse_sample_count_data = defaultdict(lambda: defaultdict(int))
    
    for file_name, (table1, table2) in prof_data_dict.items():
        damage_type = extract_damage_type(file_name)
        aligner = extract_aligner(file_name)

        if table1 is not None:
            true_data_key = f'{damage_type}{len(table1)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            rmse_table1, sample_count_table1 = compute_rmse_divergence(true_data, table1, aligner, damage_type)

            if rmse_table1 is not None:
                rmse_values_data[aligner][damage_type].append(rmse_table1)
                rmse_sample_count_data[aligner][damage_type] += sample_count_table1

        if table2 is not None:
            true_data_key = f'{damage_type}{len(table2)}.dat'
            if true_data_key in ['high5.dat', 'high3.dat', 'mid5.dat', 'mid3.dat']:
                true_data_key = "d" + true_data_key

            true_data = damage_data_dict.get(true_data_key)
            rmse_table2, sample_count_table2 = compute_rmse_divergence(true_data, table2, aligner, damage_type)

            if rmse_table2 is not None:
                rmse_values_data[aligner][damage_type].append(rmse_table2)
                rmse_sample_count_data[aligner][damage_type] += sample_count_table2

    rmse_median_data = defaultdict(lambda: defaultdict(float))
    for aligner, damage_data in rmse_values_data.items():
        for damage_type, rmse_values in damage_data.items():
            rmse_median_data[aligner][damage_type] = np.median(rmse_values)

    plot_rmse(rmse_median_data, rmse_sample_count_data, 'Median RMSE by Aligner and Damage Type', 'rmse_plot.png')
    filtered_rmse_data = filter_rmse_data_for_giraffe_and_safari(rmse_median_data)

if __name__ == "__main__":
    damage_data_path = '.'  # Your path to damage data files
    prof_data_path = 'alignments/profs'  # Your path to prof data files

    assert os.path.exists(damage_data_path), "Damage data path does not exist."
    assert os.path.exists(prof_data_path), "Prof data path does not exist."

    damage_data_files = glob.glob(os.path.join(damage_data_path, '*.dat'))
    prof_data_files = glob.glob(os.path.join(prof_data_path, '*.prof'))

    assert damage_data_files, "No damage data files found."
    assert prof_data_files, "No prof data files found."

    damage_data_dict = {os.path.basename(file_name): load_damage_data(os.path.basename(file_name)) for file_name in damage_data_files}
    prof_data_dict = {os.path.basename(file_name): load_prof_data(os.path.basename(file_name)) for file_name in prof_data_files}

    check_data(damage_data_dict, prof_data_dict)

