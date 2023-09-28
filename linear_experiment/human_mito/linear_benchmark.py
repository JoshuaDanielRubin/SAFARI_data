import csv
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def read_csv(file_path):
    data = defaultdict(lambda: defaultdict(dict))
    with open(file_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        for row in csvreader:
            damage_type, aligner_name, avg_proportion_correct, avg_proportion_incorrect = row
            data[damage_type][aligner_name] = {
                'Avg_Proportion_Correct': float(avg_proportion_correct),
                'Avg_Proportion_Incorrect': float(avg_proportion_incorrect)
            }
    return data

def plot_data(data):
    bar_width = 0.35  # Width of the bars
    for damage_type, aligner_data in data.items():
        # Rename the damage levels
        readable_damage_type = damage_type.replace("dd", "").replace("dnone", "none").replace("dsingle", "single")
        
        aligners = list(aligner_data.keys())
        avg_proportion_correct = [aligner_data[aligner]['Avg_Proportion_Correct'] for aligner in aligners]
        avg_proportion_incorrect = [aligner_data[aligner]['Avg_Proportion_Incorrect'] for aligner in aligners]
        
        safari_index = aligners.index('safari')
        
        # Create the bar plot
        plt.figure(figsize=(14, 6))
        
        # Set positions for bars
        r1 = np.arange(len(aligners))
        r2 = [x + bar_width for x in r1]
        
        # Plot both metrics
        plt.bar(r1, avg_proportion_correct, color='#56B4E9', width=bar_width, label='Proportion Correctly Mapped')
        plt.bar(r2, avg_proportion_incorrect, color='#D55E00', width=bar_width, label='Proportion Incorrectly Mapped')
        
        plt.xlabel('Aligner Name', fontsize=14)
        plt.ylabel('Proportion', fontsize=14)
        plt.title(f'Damage Level: {readable_damage_type}', fontsize=16)
        plt.xticks([r + bar_width / 2 for r in range(len(aligners))], aligners, rotation=45)
        plt.legend(fontsize=9)
        plt.grid(axis='y')
        
        plt.tight_layout()
        plt.savefig(f"Proportions_{readable_damage_type}.png")
        plt.close()

# Read data from the CSV file
file_path = 'average_proportions.csv'
data = read_csv(file_path)

# Plot the data
plot_data(data)

