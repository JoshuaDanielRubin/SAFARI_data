import os
import glob
import pandas as pd
import numpy as np
import io
import re

def clean_data(df):
    df = df.applymap(lambda x: re.search(r"([\d\.]+)", str(x)).group(1) if re.search(r"([\d\.]+)", str(x)) else np.nan)
    return df.astype(float)

def extract_damage_type(file_name):
    parts = file_name.split('_')
    for part in parts:
        if part.startswith('d') or part.startswith('dd'):
            damage_type = part.lstrip('d')
            assert damage_type in ['none', 'mid', 'high', 'single'], f'Unexpected damage type: {damage_type}'
            return damage_type
    raise ValueError(f'Unexpected file name structure, no damage type found: {file_name}')

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
        print(f'No data to check for {file_name}')
    else:
        pass #print(f'Data in {file_name} loaded successfully.')
    return data

def load_prof_data(file_name):
    file_path = os.path.join(prof_data_path, file_name)
    if os.path.getsize(file_path) == 0:
        print(f'File {file_name} is empty.')
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
        pass #print(f'Data in {file_name} loaded successfully.')
    return table1, table2

def compute_mse(true_data, estimated_data):
    if true_data is None or estimated_data is None:
        print('Missing data, cannot compute MSE.')
        return None
    try:
        true_data = clean_data(true_data)
        estimated_data = clean_data(estimated_data)

        common_columns = true_data.columns.intersection(estimated_data.columns)
        true_data = true_data[common_columns]
        estimated_data = estimated_data[common_columns]

        if true_data.empty or estimated_data.empty:
            print('Missing data, cannot compute MSE.')
            return None
        
        mse = ((true_data.values - estimated_data.values) ** 2).mean()
    except Exception as e:
        print(f'Error computing MSE: {e}')
        return None
    return mse

def check_data(damage_data_dict, prof_data_dict):
    for file_name, (table1, table2) in prof_data_dict.items():
        if table1 is None or table2 is None:
            continue
        damage_type = extract_damage_type(file_name)
        #print(f'Processing {file_name} with damage type {damage_type}')

        true_data_key = f'{damage_type}{len(table1)}.dat'
        true_data = damage_data_dict.get(true_data_key)
        
        #print(f'True Data 1:\n{true_data}')
        #print(f'Table 1:\n{table1}')
        
        mse_table1 = compute_mse(true_data, table1)
        mse_table2 = compute_mse(true_data, table2)

        if mse_table1 is not None:
            print(f'MSE for table 1: {mse_table1}')
        if mse_table2 is not None:
            print(f'MSE for table 2: {mse_table2}')

if __name__ == "__main__":
    damage_data_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito'
    prof_data_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/new_alignments/profs'

    damage_data_files = glob.glob(os.path.join(damage_data_path, '*.dat'))
    prof_data_files = glob.glob(os.path.join(prof_data_path, '*.prof'))

    damage_data_dict = {os.path.basename(file_name): load_damage_data(os.path.basename(file_name)) for file_name in damage_data_files}
    prof_data_dict = {os.path.basename(file_name): load_prof_data(os.path.basename(file_name)) for file_name in prof_data_files}

    check_data(damage_data_dict, prof_data_dict)

