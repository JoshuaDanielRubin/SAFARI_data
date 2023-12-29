
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def calculate_precision_recall(df):
    precision = df['TP'] / (df['TP'] + df['FP'])
    recall = df['TP'] / (df['TP'] + df['FN'])
    return precision, recall

def main(file_path):
    # Load data
    data = pd.read_csv(file_path)

    # Filter data for 'giraffe' and 'SAFARI'
    giraffe_data = data[data['Aligner_Name'].str.lower() == 'giraffe']
    safari_data = data[data['Aligner_Name'].str.upper() == 'SAFARI']

    # Damage types
    damage_types = data['Damage_Type'].unique()

    # Colorblind-friendly colors
    colors = ['orange', 'teal', 'magenta', 'lime']

    # Line styles
    line_styles = {'giraffe': '-', 'SAFARI': '--'}

    # Line thickness
    line_thickness = 2.5

    # Prepare plot
    plt.figure(figsize=(14, 10))

    # Plotting
    for i, damage_type in enumerate(damage_types):
        # Filter data by damage type
        giraffe_dt_data = giraffe_data[giraffe_data['Damage_Type'] == damage_type]
        safari_dt_data = safari_data[safari_data['Damage_Type'] == damage_type]

        # Calculate precision and recall
        giraffe_precision, giraffe_recall = calculate_precision_recall(giraffe_dt_data)
        safari_precision, safari_recall = calculate_precision_recall(safari_dt_data)

        # Fit and plot curves for giraffe data
        if len(giraffe_precision) > 1 and len(giraffe_recall) > 1:
            giraffe_curve_fit = np.polyfit(giraffe_recall, giraffe_precision, 2)
            giraffe_fit_fn = np.poly1d(giraffe_curve_fit)
            recall_range = np.linspace(min(giraffe_recall), max(giraffe_recall), 100)
            plt.plot(recall_range, giraffe_fit_fn(recall_range), line_styles['giraffe'], color=colors[i], linewidth=line_thickness, label=f'giraffe - {damage_type}')

        # Fit and plot curves for SAFARI data
        if len(safari_precision) > 1 and len(safari_recall) > 1:
            safari_curve_fit = np.polyfit(safari_recall, safari_precision, 2)
            safari_fit_fn = np.poly1d(safari_curve_fit)
            recall_range = np.linspace(min(safari_recall), max(safari_recall), 100)
            plt.plot(recall_range, safari_fit_fn(recall_range), line_styles['SAFARI'], color=colors[i], linewidth=line_thickness, label=f'SAFARI - {damage_type}')

    # Adding labels and title
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Tradeoff: giraffe vs SAFARI (Stratified by Damage Type)', fontsize=22)
    plt.legend(loc='lower left')
    plt.grid(True)

    # Save the plot with high resolution
    plt.savefig('precision_recall_tradeoff.png', dpi=300)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python this_script.py <path_to_alignment_stats.csv>")
    else:
        main(sys.argv[1])

