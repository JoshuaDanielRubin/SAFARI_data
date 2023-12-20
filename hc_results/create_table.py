import os
import pandas as pd
import re  # Added for extracting replicate numbers

# Specify the directory containing your log files
log_dir = "."

def extract_info_from_filename(filename):
    """Extract sample name, subsampling rate, correction status, and replicate number from the filename."""
    
    # Extracting sample name using regex
    sample_name_match = re.search(r'^(.*?)_', filename)
    sample_name = sample_name_match.group(1).split(".")[0] if sample_name_match else filename.split('.')[0]
    
    subsampling_rate = next((segment.split('x')[0] for segment in filename.split('_') if 'x' in segment), None)
    correction_status = filename.split('.')[-2]
    
    # Extracting replicate number using regex
    replicate_match = re.search(r'_replicate_(\d+)', filename)
    replicate_number = replicate_match.group(1) if replicate_match else "1"
    
    return sample_name, float(subsampling_rate), correction_status, replicate_number

def extract_info_from_file(filepath):
    """Extract the haplogroup and number of reads from the file content."""
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('stdin'):
                parts = line.strip().split()
                haplogroup = parts[1]
                reads = int(parts[2])
                return haplogroup, reads
    return None, None

# Create a dictionary mapping sample names to combined predictions
full_coverage_predictions = {
    "DA100": "C4b1/C4b1",
    "DA101": "U5a1b1e/U5a1b1e",
      "DA171": "H2a1/H2a1",
    "DA15": "C4d/C4d",
    "I10899": "U5b/U5b",
    "I11300": "J2a1a1/J2a1a1",
    "I8132": "D1a2/H2a2a1",
    "I8569": "H1ah/H1ah",
    "I7645": "R7b2/H2a2a1",
    "I7646": "H1e1c/H2a2a1",
    "STR393b": "H5a1/H5a1",
    "NW54": "C4a1a1/C4a1a1",
    "AED204": "X2b11/X2b+226",
    "Alh10": "I1/I1",
    "FN2": "H3/H3",
    "AED92b": "U4a1/U4a1",
    "STR266b": "J1c/J1c5",
    "Vim2b": "H7/H7",
    "STR491": "T2b/T2b",
    "STR486": "T2b/T2b",
}

data = []

# Iterate through all log files in the specified directory
for filename in os.listdir(log_dir):
    if filename.endswith('.log'):
        sample_name, subsampling_rate, correction_status, replicate_number = extract_info_from_filename(filename)
        filepath = os.path.join(log_dir, filename)
        haplogroup, reads = extract_info_from_file(filepath)

        data.append({
            'Sample Name': sample_name,
            'Subsampling Rate': subsampling_rate,
            'Correction Status': correction_status,
            'Haplogroup': haplogroup,
            'Reads': reads,
            'Replicate': replicate_number,
            'Full Coverage Prediction': full_coverage_predictions.get(sample_name, "N/A")
        })

# Create a DataFrame from the data
df = pd.DataFrame(data)

# Pivot the table to have separate columns for corrected and uncorrected data
pivot_df = df.pivot_table(index=['Sample Name', 'Subsampling Rate', 'Full Coverage Prediction', 'Replicate'], columns='Correction Status',
                          values=['Haplogroup', 'Reads'], aggfunc='first')

# Flatten the MultiIndex to have single-level columns
pivot_df.columns = ['_'.join(col).strip() for col in pivot_df.columns.values]

# Reset the index to bring Sample Name, Subsampling Rate, and Replicate back as columns
pivot_df.reset_index(inplace=True)

# Adding 'X' to Subsampling Rate values
pivot_df['Subsampling Rate'] = pivot_df['Subsampling Rate'].astype(str) + 'X'

# Renaming columns to shorter names
pivot_df.rename(columns={
    'Haplogroup_corrected': 'HG_corrected',
    'Haplogroup_uncorrected': 'HG_uncorrected',
    #'Reads_corrected': '# Reads_corrected',
    #'Reads_uncorrected': '# Reads_uncorrected',
    'Subsampling Rate': 'Rate',
    'Full Coverage Prediction': 'Full_Coverage_Prediction'
}, inplace=True)

# Sort the DataFrame by subsampling rate
pivot_df = pivot_df.sort_values(by=['Rate', 'Sample Name', 'Replicate'])

# Save the DataFrame to a LaTeX file
with open('table.txt', 'w') as f:
    latex_string = pivot_df.to_latex(index=False, na_rep='N/A', longtable=True)
    latex_string += '\n\\caption{This is the caption for the table.}'
    f.write(latex_string)

print("Data aggregation complete. The data has been saved to 'table.txt'.")

