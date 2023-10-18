import pandas as pd
import os

# Directory containing the TSV files
directory = "/home/projects/MAAG/Magpie/Magpie/linear_experiment/human_mito/benchmarks"

# Lists to store dataframes for each category
safari_dfs = []
giraffe_dfs = []

# List of critical metrics
important_metrics = ['s', 'max_rss', 'mean_load']

# Iterate over each TSV file and read into a DataFrame
for filename in os.listdir(directory):
    path = os.path.join(directory, filename)
    if filename.endswith(".tsv"):
        df = pd.read_csv(path, sep='\t')[important_metrics]  # Filter by important metrics
        if "safari" in filename:
            safari_dfs.append(df)
        elif "giraffe" in filename:
            giraffe_dfs.append(df)

# Aggregate the dataframes for each category (assuming you want mean values)
safari_aggregated = pd.concat(safari_dfs, ignore_index=True).mean()
giraffe_aggregated = pd.concat(giraffe_dfs, ignore_index=True).mean()

# Prepare a DataFrame for LaTeX conversion
combined_df = pd.DataFrame({
    "Metrics": safari_aggregated.index,
    "Safari": safari_aggregated.values,
    "Giraffe": giraffe_aggregated.values
})

# Convert to LaTeX with a suitable caption and save to a file
latex_code = combined_df.to_latex(index=False, 
                                  caption="Comparison of key performance metrics between Safari and Giraffe.", 
                                  label="tab:key_metrics_comparison")
with open(os.path.join(directory, "key_metrics_comparison_table.tex"), "w") as f:
    f.write(latex_code)

