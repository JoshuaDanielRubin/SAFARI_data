import os
import glob
import pandas as pd
import numpy as np
import io
import re
import matplotlib.pyplot as plt
from collections import defaultdict

def clean_data(df):
    df = df.applymap(lambda x: re.search(r"([\d\.]+)", str(x)).group(1) if re.search(r"([\d\.]+)", str(x)) else np.nan)
    return df.astype(float)

def extract_damage_type(file_name):
    parts = file_name.split('_')
    for part in parts:
        if part.startswith('d') or part.startswith('dd'):
            damage_type = part.lstrip('d').lstrip('d')  # Adjusted this line
            assert damage_type in ['none', 'mid', 'high', 'single'], f'Unexpected damage type: {damage_type}'
            return damage_type
    raise ValueError(f'Unexpected file name structure, no damage type found: {file_name}')

def extract_aligner(file_name):
    return file_name.split('_')[-1].split('.')[0]

def load_damage_data(file_name):
    file_path = os.path.join(damage_data_path, file_name)
    if os.path.getsize(file_path) == 0:
        print(f'File {file_name} is empty.')
        return None
    try:
        data = pd.read_csv(file_path, delimiter='\t', index_col=0)
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None
    if data.empty:
        print(f'No data in file {file_name}.')
    else:
        pass
    return data

def load_prof_data(file_name):
    file_path = os.path.join(prof_data_path, file_name)
    if os.path.getsize(file_path) == 0:
        #print(f'File {file_name} is empty.')
        return None, None
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            second_header_index = len(lines) // 2
            table1_lines = lines[:second_header_index]
            table2_lines = lines[second_header_index:]

        table1 = pd.read_csv(io.StringIO(''.join(table1_lines)), delimiter='\t', header=0)
        table2 = pd.read_csv(io.StringIO(''.join(table2_lines)), delimiter='\t', header=0)
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None, None

    if table1.empty or table2.empty:
        print(f'No data to check for {file_name}')
    else:
        pass
    return table1, table2

def compute_mse(true_data, estimated_data):
    if true_data is None:
        print('True data is missing.')
        return None
    if estimated_data is None:
        print('Estimated data is missing.')
        return None
    try:
        true_data = clean_data(true_data)
        estimated_data = clean_data(estimated_data)

        common_columns = true_data.columns.intersection(estimated_data.columns)
        true_data = true_data[common_columns]
        estimated_data = estimated_data[common_columns]

        if true_data.empty:
            print('True data is empty.')
            return None
        if estimated_data.empty:
            print('Estimated data is empty.')
            return None
        
        mse = ((true_data.values - estimated_data.values) ** 2).mean()
    except Exception as e:
        print(f'Error computing MSE: {e}')
        return None
    return mse

def filter_mse_data_for_giraffe_and_safari(mse_data):
    filtered_mse_data = defaultdict(lambda: defaultdict(float))
    for aligner in ['giraffe', 'safari']:
        if aligner in mse_data:
            filtered_mse_data[aligner] = mse_data[aligner]
    return filtered_mse_data

def plot_mse(mse_data, plot_title, save_file_name):
    colorblind_colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', '#CA9161', '#FBAFE4']
    aligners = list(mse_data.keys())
    damage_types = list(mse_data[aligners[0]].keys())
    
    x = np.arange(len(aligners))
    width = 0.2

    fig, ax = plt.subplots()
    for i, damage_type in enumerate(damage_types):
        mse_values = [mse_data[aligner][damage_type] for aligner in aligners]
        ax.bar(x + i*width, mse_values, width, label=damage_type, color=colorblind_colors[i])
    
    ax.set_xlabel('Aligner')
    ax.set_ylabel('MSE')
    ax.set_title(plot_title)
    ax.set_xticks(x + width*(len(damage_types)-1)/2)
    ax.set_xticklabels(aligners)
    ax.legend()

    fig.tight_layout()
    plt.savefig(save_file_name)
    
    # Print MSE values
    for aligner in aligners:
        for damage_type in damage_types:
            print(f'MSE for {aligner} with damage type {damage_type}: {mse_data[aligner][damage_type]}')

def check_data(damage_data_dict, prof_data_dict):
    mse_data = defaultdict(lambda: defaultdict(float))
    for file_name, (table1, table2) in prof_data_dict.items():
        if table1 is None or table2 is None:
            continue
        damage_type = extract_damage_type(file_name)
        aligner = extract_aligner(file_name)

        true_data_key = f'{damage_type}{len(table1)}.dat'

        if true_data_key == 'high5.dat':
            true_data_key= "dhigh5.dat"
        if true_data_key == 'high3.dat':
            true_data_key= "dhigh3.dat"
        if true_data_key == 'mid5.dat':
            true_data_key= "dmid5.dat"
        if true_data_key == 'mid3.dat':
            true_data_key= "dmid3.dat"

        true_data = damage_data_dict.get(true_data_key)

        mse_table1 = compute_mse(true_data, table1)
        mse_table2 = compute_mse(true_data, table2)

        if mse_table1 is not None:
            mse_data[aligner][damage_type] = mse_table1
        if mse_table2 is not None:
            mse_data[aligner][damage_type] = mse_table2

    plot_mse(mse_data, 'MSE by aligner and damage type', 'mse_plot.png')
    filtered_mse_data = filter_mse_data_for_giraffe_and_safari(mse_data)
    plot_mse(filtered_mse_data, 'MSE comparison between Giraffe and Safari', 'mse_plot_giraffe_safari.png')

if __name__ == "__main__":
    damage_data_path = '.'
    prof_data_path = 'new_alignments/profs'

    damage_data_files = glob.glob(os.path.join(damage_data_path, '*.dat'))
    prof_data_files = glob.glob(os.path.join(prof_data_path, '*.prof'))

    damage_data_dict = {os.path.basename(file_name): load_damage_data(os.path.basename(file_name)) for file_name in damage_data_files}
    prof_data_dict = {os.path.basename(file_name): load_prof_data(os.path.basename(file_name)) for file_name in prof_data_files}

    check_data(damage_data_dict, prof_data_dict)

