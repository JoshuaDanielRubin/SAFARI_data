import pandas as pd
from io import StringIO

pd.set_option('display.max_rows', None)

def process_latex_data_for_pandas(latex_table):
    lines = latex_table.split('\n')
    data_lines = [line for line in lines if line.strip() and not line.strip().startswith('\\')]
    data_lines = [line for line in data_lines if not line.strip().startswith('end')]
    processed_lines = []
    for line in data_lines:
        cells = [cell.replace('\\', '').replace('BETTER', '').replace('WORSE', '').strip() for cell in line.split('&') if cell.strip()]
        processed_line = ','.join(cells)
        processed_lines.append(processed_line)
    csv_table = '\n'.join(processed_lines)
    return csv_table

# Replace 'file_path' with the path to your LaTeX table file
file_path = 'table.txt'

with open(file_path, 'r') as file:
    full_latex_table = file.read()

csv_table_improved = process_latex_data_for_pandas(full_latex_table)
df_full_improved = pd.read_csv(StringIO(csv_table_improved))
# Convert to numeric, coerce non-numeric to NaN
df_full_improved['Reads_corrected'] = pd.to_numeric(df_full_improved['Reads_corrected'], errors='coerce')
df_full_improved['Reads_uncorrected'] = pd.to_numeric(df_full_improved['Reads_uncorrected'], errors='coerce')

# Assuming 'Reads_corrected' is the 7th column and 'Reads_uncorrected' is the 8th column
# Adjust the indices as needed based on your actual data
mean_corrected = df_full_improved.iloc[1:, 6].dropna().mean()  # Start from the second row
mean_uncorrected = df_full_improved.iloc[1:, 7].dropna().mean()  # Start from the second row

print(f"mean Reads Corrected: {mean_corrected}")
print(f"mean Reads Uncorrected: {mean_uncorrected}")

# Further processing and analysis can be done on df_full_improved as needed

