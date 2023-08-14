import pandas as pd
import matplotlib.pyplot as plt
import os

# Colorblind-friendly colors
colors_minimizer = ['#A6CEE3', '#1F78B4']  # Light Blue, Dark Blue
colors_rymer = ['#FB9A99', '#E31A1C']  # Light Red, Dark Red

# Function to generate stratified plots by a given column
def generate_stratified_plot_by_column(column_name, save_name):
    unique_values = df[column_name].unique()
    num_values = len(unique_values)
    
    fig, axes = plt.subplots(nrows=num_values, figsize=(14, 6 * num_values), sharex=True, sharey=True)

    # Ensure axes is always a list for consistent indexing
    if num_values == 1:
        axes = [axes]

    for i, value in enumerate(sorted(unique_values)):
        ax = axes[i]
        df_filtered = df[df[column_name] == value]
        
        df_grouped = df_filtered.groupby(['N', 'Minimizer sketch']).mean()[['Total minimizer seeds', 'Total rymer seeds']]
        
        # Plotting for each unique value
        df_grouped.unstack()['Total minimizer seeds'].plot(kind='bar', position=0, width=0.4, ax=ax, color=colors_minimizer)
        df_grouped.unstack()['Total rymer seeds'].plot(kind='bar', position=1, width=0.4, ax=ax, color=colors_rymer)

        # Add title, labels, and other plot properties for each subplot
        ax.set_title(f'Distribution of Minimizer and Rymer Seed Counts by Sketch Type ({column_name}={value})')
        ax.set_xlabel('Number of Reads')
        ax.set_ylabel('Seed Counts')
        ax.set_xticklabels(df_grouped.index.levels[0], rotation=45)
        ax.legend(['Minimizer (Naive)', 'Minimizer (Unique)', 'Rymer (Naive)', 'Rymer (Unique)'])
        ax.grid(axis='y')
    
    plt.tight_layout()
    plt.savefig(f'plots/{save_name}.png')

# Create a directory for plots if it doesn't exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Read the data
df = pd.read_csv('results.tsv', sep='\t')

# Call the functions to generate and save plots
# Uncomment the following lines to generate the plots:
generate_stratified_plot_by_column('delta', 'stratified_plot_by_delta')
generate_stratified_plot_by_column('W', 'stratified_plot_by_W')
generate_stratified_plot_by_column('L', 'stratified_plot_by_read_length')

