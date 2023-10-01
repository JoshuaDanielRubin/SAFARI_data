import os
import glob
import pandas as pd
import numpy as np

def load_damage_data(file_name):
    file_path = os.path.join(damage_data_path, file_name)
    if os.path.getsize(file_path) == 0:
        print(f'File {file_name} is empty.')
        return None
    try:
        data = pd.read_csv(file_path, header=None)
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None
    if data.empty:
        print(f'No data to check for {file_name}')
    else:
        print(f'Data in {file_name} loaded successfully.')
    return data

def load_prof_data(file_name):
    file_path = os.path.join(prof_data_path, file_name)
    if os.path.getsize(file_path) == 0:
        print(f'File {file_name} is empty.')
        return None, None
    try:
        table1 = pd.read_csv(file_path, header=None, skiprows=1, nrows=5)
        table2 = pd.read_csv(file_path, header=None, skiprows=7)
    except Exception as e:
        print(f'Error reading {file_name}: {e}')
        return None, None
    if table1.empty or table2.empty:
        print(f'No data to check for {file_name}')
    else:
        print(f'Data in {file_name} loaded successfully.')
    return table1, table2

def calculate_accuracy(predicted_data, ground_truth_data):
    if predicted_data is None or ground_truth_data is None or predicted_data.shape != ground_truth_data.shape:
        return None
    sum_of_absolute_differences = np.sum(np.abs(predicted_data.values - ground_truth_data.values))
    total_elements = np.prod(predicted_data.shape)
    accuracy = 1 - (sum_of_absolute_differences / total_elements)
    return accuracy

def check_data(damage_data, prof_data):
    accuracy_results = {}

    for damage_file, damage_df in damage_data.items():
        prof_file = damage_file.replace('.dat', '.prof')
        prof_tables = prof_data.get(prof_file)

        if prof_tables is None:
            print(f"No prof data found for {damage_file}")
            continue

        table1, table2 = prof_tables
        accuracy_results[(damage_file, 'single3.dat')] = calculate_accuracy(table1, damage_df)
        accuracy_results[(damage_file, 'single5.dat')] = calculate_accuracy(table2, damage_df)

    for key, accuracy in accuracy_results.items():
        print(f'Accuracy for {key[0]} with {key[1]}: {accuracy}')

if __name__ == "__main__":
    damage_data_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito'
    prof_data_path = '/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/new_alignments/profs'

    damage_data_files = glob.glob(os.path.join(damage_data_path, '*.dat'))
    prof_data_files = glob.glob(os.path.join(prof_data_path, '*.prof'))

    damage_data = {os.path.basename(file_name): load_damage_data(os.path.basename(file_name)) for file_name in damage_data_files}
    prof_data = {os.path.basename(file_name): load_prof_data(os.path.basename(file_name)) for file_name in prof_data_files}

    check_data(damage_data, prof_data)

