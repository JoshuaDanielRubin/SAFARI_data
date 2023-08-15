import pandas as pd
import matplotlib.pyplot as plt
import os

# Ensure the 'plots' directory exists
if not os.path.exists('plots'):
    os.makedirs('plots')

# Load the data from the TSV file
data = pd.read_csv("results.tsv", sep="\t")

# Calculate the expected Rymer rescue rate
data['Expected Rate'] = (data['delta']**2) / (1 + 2*data['delta'] + data['delta']**2)

# Group by delta and compute the mean and standard deviation of the Correct Rescue Rate
grouped_data = data.groupby('delta').agg({
    'Correct Rescue Rate': ['mean', 'std'],
    'Expected Rate': 'first'  # Expected rate is constant for each delta value, so taking the first value is sufficient
}).reset_index()

# Plotting with error bars and colorblind-friendly colors
plt.figure(figsize=(10, 6))
plt.errorbar(grouped_data['delta'], grouped_data[('Correct Rescue Rate', 'mean')],
             yerr=grouped_data[('Correct Rescue Rate', 'std')],
             fmt='o', color='#0173B2', label='Actual Rate (Mean ± SD)', capsize=5)
plt.scatter(grouped_data['delta'], grouped_data['Expected Rate'], color='#DE8F05', label='Expected Rate', s=80)
plt.title('Comparison of Actual vs Expected Rymer Rescue Rate (Grouped by Delta)')
plt.xlabel('Deamination Rate')
plt.ylabel('Correct Rescue Rate')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the plot to the specified directory
plt.savefig("plots/fit.png")

