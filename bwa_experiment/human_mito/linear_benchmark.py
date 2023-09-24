
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Data
data = [(0.01, 'safari', 1.2, 20.76), (0.01, 'bwa-aln_anc', 1.2, 9.11), (0.01, 'bwa-aln', 1.2, 6.09), (0.01, 'bwa-mem', 1.08, 7.03), (0.01, 'BBMap', 1.2, 14.95), (0.01, 'shrimp', 1.2, 23.71), (0.01, 'bowtie2', 1.2, 5.39), (0.03, 'safari', 2.81, 18.06), (0.03, 'bwa-aln_anc', 2.79, 7.88), (0.03, 'bwa-aln', 2.79, 5.98), (0.03, 'bwa-mem', 2.51, 6.85), (0.03, 'BBMap', 2.8, 15.82), (0.03, 'shrimp', 2.8, 22.48), (0.03, 'bowtie2', 2.78, 4.37), (0.05, 'safari', 4.1, 17.2), (0.05, 'bwa-aln_anc', 4.05, 7.22), (0.05, 'bwa-aln', 4.05, 5.57), (0.05, 'bwa-mem', 3.65, 6.23), (0.05, 'BBMap', 4.09, 14.95), (0.05, 'shrimp', 4.09, 21.23), (0.05, 'bowtie2', 4.04, 4.03), (0.07, 'safari', 6.59, 17.89), (0.07, 'bwa-aln_anc', 6.51, 7.57), (0.07, 'bwa-aln', 6.5, 5.81), (0.07, 'bwa-mem', 5.86, 6.42), (0.07, 'BBMap', 6.57, 15.12), (0.07, 'shrimp', 6.57, 21.44), (0.07, 'bowtie2', 6.47, 3.97), (0.1, 'safari', 9.48, 17.78), (0.1, 'bwa-aln_anc', 9.39, 6.66), (0.1, 'bwa-aln', 9.37, 4.72), (0.1, 'bwa-mem', 8.44, 5.97), (0.1, 'BBMap', 9.46, 14.64), (0.1, 'shrimp', 9.46, 21.65), (0.1, 'bowtie2', 9.35, 3.35), (0.3, 'safari', 32.6, 17.82), (0.3, 'bwa-aln_anc', 32.32, 7.27), (0.3, 'bwa-aln', 32.18, 4.99), (0.3, 'bwa-mem', 28.78, 6.1), (0.3, 'BBMap', 32.51, 15.2), (0.3, 'shrimp', 32.51, 21.93), (0.3, 'bowtie2', 32.12, 3.3), (0.5, 'safari', 51.49, 17.65), (0.5, 'bwa-aln_anc', 51.15, 7.56), (0.5, 'bwa-aln', 50.93, 5.01), (0.5, 'bwa-mem', 45.52, 6.18), (0.5, 'BBMap', 51.37, 14.83), (0.5, 'shrimp', 51.37, 22.09), (0.5, 'bowtie2', 50.66, 3.34), (0.7, 'safari', 70.48, 17.91), (0.7, 'bwa-aln_anc', 69.99, 7.69), (0.7, 'bwa-aln', 69.65, 5.36), (0.7, 'bwa-mem', 62.28, 6.44), (0.7, 'BBMap', 70.3, 15.07), (0.7, 'shrimp', 70.29, 22.02), (0.7, 'bowtie2', 69.29, 3.61), (0.9, 'safari', 90.03, 17.99), (0.9, 'bwa-aln_anc', 89.36, 7.76), (0.9, 'bwa-aln', 88.96, 5.4), (0.9, 'bwa-mem', 79.55, 6.51), (0.9, 'BBMap', 89.78, 15.02), (0.9, 'shrimp', 89.77, 22.03), (0.9, 'bowtie2', 88.5, 3.77)]
columns = ['Sampling Rate', 'Aligner', 'Correct (%)', 'Incorrect (%)']

# DataFrame
df = pd.DataFrame(data, columns=columns)
df.set_index(['Sampling Rate', 'Aligner'], inplace=True)

# Function to calculate percentage difference
def calculate_percentage_diff_positive(group):
    safari_values = group.loc[group.index.get_level_values('Aligner') == 'safari']
    return safari_values.values - group  # Safari values - Tool values to get positive differences

# Group by Sampling Rate and apply the function
df_diff_positive = df.groupby('Sampling Rate').apply(calculate_percentage_diff_positive)
df_diff_positive = df_diff_positive[df_diff_positive.index.get_level_values('Aligner') != 'safari'].reset_index()

# Set the color palette to a color-blind friendly palette
sns.set_palette("colorblind")
sns.set(font_scale=0.9)  # Set font size

# Get the unique sampling rates
unique_sampling_rates = df_diff_positive['Sampling Rate'].unique()

# Determine the number of rows required for a 2-column layout
num_rows = (len(unique_sampling_rates) + 1) // 2

# Create subplots with 2 columns
fig, axs = plt.subplots(num_rows, 2, figsize=(15, 5 * num_rows))
axs = axs.flatten()  # Flatten the axes array for easier iteration

# Iterate through each sampling rate and create a bar plot for Correct (%)
for i, rate in enumerate(unique_sampling_rates):
    df_subset = df_diff_positive[df_diff_positive['Sampling Rate'] == rate]
    sns.barplot(data=df_subset, x='Aligner', y='Correct (%)', ax=axs[i])

    # Set titles and labels with multiple lines for clarity
    title = (
        f'Sampling Rate {rate}:\n'
        'Percentage Difference in Number of Mappable Reads Correctly Mapped\n'
        '(Comparing Safari to Other Tools)'
    )
    axs[i].set_title(title.format(rate=rate))
    axs[i].set_ylabel('Percentage Difference')
    axs[i].tick_params(axis='x', rotation=45)

# Hide any unused subplots
for i in range(len(unique_sampling_rates), len(axs)):
    axs[i].axis('off')

# Adjust the layout and increase the spacing between subplots
plt.tight_layout()
plt.subplots_adjust(wspace=0.3)  # Adjusting the space between columns

# Show the plot
plt.savefig("linear_benchmark.png")
plt.close()

