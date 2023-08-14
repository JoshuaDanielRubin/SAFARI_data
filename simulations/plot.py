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
generate_stratified_plot_by_column('delta', 'stratified_plot_by_delta')
generate_stratified_plot_by_column('W', 'stratified_plot_by_W')
generate_stratified_plot_by_column('L', 'stratified_plot_by_read_length')

def plot_side_by_side_subplots(df):
    """Plot two side-by-side subplots comparing precision and correct rescue rate for both rymers and minimizers."""
    # Set up the figure and axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Define the bar width and positions
    bar_width = 0.35
    indices = range(len(df))
    
    # First subplot for precision
    ax1.bar(indices, df["Minimizer precision"], bar_width, label='Minimizer Precision', color=colors_minimizer[1])
    ax1.bar([i + bar_width for i in indices], df["Rymer precision"], bar_width, label='Rymer Precision', color=colors_rymer[1])
    ax1.set_title('Precision Comparison')
    ax1.set_xticks([i + bar_width/2 for i in indices])
    ax1.set_xticklabels(indices)  # Using row indices for labeling
    ax1.legend()
    ax1.set_xlabel('Data Points')
    ax1.set_ylabel('Precision')
    
    # Second subplot for correct rescue rate
    ax2.bar(indices, df["Correct Rescue Rate"], bar_width, label='Minimizer', color=colors_minimizer[0])
    ax2.bar([i + bar_width for i in indices], df["Correct Rescue Rate"], bar_width, label='Rymer', color=colors_rymer[0])
    ax2.set_title('Correct Rescue Rate Comparison')
    ax2.set_xticks([i + bar_width/2 for i in indices])
    ax2.set_xticklabels(indices)  # Using row indices for labeling
    ax2.legend()
    ax2.set_xlabel('Data Points')
    ax2.set_ylabel('Correct Rescue Rate')
    
    # Display the plots
    plt.tight_layout()
    plt.savefig('plots/side_by_side_comparison.png')
    plt.show()

# Call the new function after reading the data
plot_side_by_side_subplots(df)

