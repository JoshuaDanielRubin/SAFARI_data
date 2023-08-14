import pandas as pd
import matplotlib.pyplot as plt

def generate_plot():
    # 1. Read the data from results.tsv
    df = pd.read_csv('results.tsv', sep='\t')

    # 2. Generate the grouped bar plot
    plt.figure(figsize=(12, 8))
    df_grouped = df.groupby(['N', 'Minimizer sketch']).mean()[['Total minimizer seeds', 'Total rymer seeds']].unstack()
    df_grouped['Total minimizer seeds'].plot(kind='bar', position=0, width=0.4, ax=plt.gca(), color=['#FF9999', '#66B2FF'])
    df_grouped['Total rymer seeds'].plot(kind='bar', position=1, width=0.4, ax=plt.gca(), color=['#FF6666', '#3399FF'])

    # 3. Add title, labels, and other plot properties
    plt.title('Distribution of Minimizer and Rymer Seed Counts by Sketch Type')
    plt.xlabel('Value of N')
    plt.ylabel('Seed Counts')
    plt.xticks(rotation=45)
    plt.legend(['Minimizer (Naive)', 'Minimizer (Unique)', 'Rymer (Naive)', 'Rymer (Unique)'])
    plt.tight_layout()
    plt.grid(axis='y')

    # 4. Save the plot to seed_counts_distribution.png
    plt.savefig('seed_counts_distribution.png')

# Call the function
generate_plot()

