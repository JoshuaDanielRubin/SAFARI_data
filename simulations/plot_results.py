import pandas as pd
import matplotlib.pyplot as plt

# Load the TSV file into a DataFrame
df = pd.read_csv('results.tsv', sep='\t')

# Define the combinations of variables to plot
combinations = [
    ('Rymer spuriousness', 'Rymer recovery rate'),
]

# Create a scatter plot for each combination
for x, y in combinations:
    plt.figure(figsize=(10, 6))
    plt.scatter(df[x], df[y])
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f'{y} vs {x}')
    plt.grid(True)

    # Save the plot to a file
    plt.savefig(f'{x}_vs_{y}.png')

plt.show()

