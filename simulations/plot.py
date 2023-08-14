import pandas as pd
import matplotlib.pyplot as plt
import os

# Colorblind-friendly colors
colors_minimizer = ['#A6CEE3', '#1F78B4']  # Light Blue, Dark Blue
colors_rymer = ['#FB9A99', '#E31A1C']  # Light Red, Dark Red

def generate_single_plot_save(k):
    df_k = df[df['k'] == k]
    df_grouped = df_k.groupby(['N', 'Minimizer sketch']).mean()[['Total minimizer seeds', 'Total rymer seeds']]
    
    plt.figure(figsize=(12, 6))
    df_grouped.unstack()['Total minimizer seeds'].plot(kind='bar', position=0, width=0.4, color=colors_minimizer)
    df_grouped.unstack()['Total rymer seeds'].plot(kind='bar', position=1, width=0.4, color=colors_rymer)

    # Add title, labels, and other plot properties
    plt.title(f'Distribution of Minimizer and Rymer Seed Counts by Sketch Type (k={k})')
    plt.xlabel('Number of Reads')
    plt.ylabel('Seed Counts')
    plt.xticks(rotation=45)
    plt.legend(['Minimizer (Naive)', 'Minimizer (Unique)', 'Rymer (Naive)', 'Rymer (Unique)'])
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f'plots/single_plot_k_{k}.png')

def generate_stratified_plot_by_k_save():
    unique_kmer_sizes = df['k'].unique()
    num_kmers = len(unique_kmer_sizes)
    
    fig, axes = plt.subplots(nrows=num_kmers, figsize=(14, 6 * num_kmers), sharex=True, sharey=True)

    # Ensure axes is always a list for consistent indexing
    if num_kmers == 1:
        axes = [axes]

    for i, k in enumerate(sorted(unique_kmer_sizes)):
        ax = axes[i]
        df_k = df[df['k'] == k]
        
        df_grouped = df_k.groupby(['N', 'Minimizer sketch']).mean()[['Total minimizer seeds', 'Total rymer seeds']]
        
        # Plotting for each k-mer size
        df_grouped.unstack()['Total minimizer seeds'].plot(kind='bar', position=0, width=0.4, ax=ax, color=colors_minimizer)
        df_grouped.unstack()['Total rymer seeds'].plot(kind='bar', position=1, width=0.4, ax=ax, color=colors_rymer)

        # Add title, labels, and other plot properties for each subplot
        ax.set_title(f'Distribution of Minimizer and Rymer Seed Counts by Sketch Type (k={k})')
        ax.set_xlabel('Number of Reads')
        ax.set_ylabel('Seed Counts')
        ax.set_xticklabels(df_grouped.index.levels[0], rotation=45)
        ax.legend(['Minimizer (Naive)', 'Minimizer (Unique)', 'Rymer (Naive)', 'Rymer (Unique)'])
        ax.grid(axis='y')
    
    plt.tight_layout()
    plt.savefig('plots/stratified_plot_by_k.png')

# Create a directory for plots if it doesn't exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Read the data
df = pd.read_csv('results.tsv', sep='\t')

# Call the functions to generate and save plots
# Uncomment the following lines to generate the plots:
#generate_single_plot_save(<specific k value>)
generate_stratified_plot_by_k_save()

# Replace <specific k value> with a k-mer size if you want to generate a single plot for that k-mer size.

