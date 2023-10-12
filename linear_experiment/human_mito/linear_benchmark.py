import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tabulate import tabulate

def process_data(df):
    required_columns = ['Damage_Type', 'Aligner_Name', 'Mapped_to_MT', 'Mapped_to_MT_Correct_Location', 'Mapped_to_MT_Correct_Location_MQ>30', 'Unmapped_Reads']
    assert set(required_columns).issubset(df.columns), "Not all required columns are present in the dataframe."

    custom_order_df = df.copy()

    custom_order_df['Mapped_NOT_Correctly'] = custom_order_df['Mapped_to_MT'] - custom_order_df['Mapped_to_MT_Correct_Location']

    damage_type_order = ['none', 'dmid', 'dhigh', 'single']
    damage_type_rename = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single'}

    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].astype('category')
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].cat.reorder_categories(damage_type_order, ordered=True)
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].map(damage_type_rename)

    aligner_rename = {
        'safari': 'SAFARI',
        'mem': 'BWA-MEM',
        'aln_anc': 'BWA-aln (anc)',
        'aln': 'BWA-aln',
        'bb': 'BBMap',
        'shrimp': 'SHRiMP'
    }
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].replace(aligner_rename)

    custom_aligner_order = ['SAFARI', 'giraffe', 'BWA-MEM', 'BWA-aln (anc)', 'BWA-aln', 'BBMap', 'SHRiMP']
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].astype('category')
    custom_order_df['Aligner_Name'].cat.set_categories(custom_aligner_order, ordered=True, inplace=True)
    custom_order_df.sort_values('Aligner_Name', inplace=True)

     # Calculate the medians
    metrics = [
        "Mapped_to_MT_Correct_Location",
        "Mapped_NOT_Correctly",
        "Mapped_to_MT_Correct_Location_MQ>30",
        "Unmapped_Reads"
    ]
    
    # Instead of using the transform, let's explicitly compute the median and then merge it back to the dataframe
    medians_df = custom_order_df.groupby(['Aligner_Name', 'Damage_Type'])[metrics].median().reset_index()
    
    custom_order_df = custom_order_df.drop(columns=metrics).drop_duplicates(subset=['Aligner_Name', 'Damage_Type'])
    custom_order_df = pd.merge(custom_order_df, medians_df, on=['Aligner_Name', 'Damage_Type'], how='left')
    
    return custom_order_df


def create_new_plot(df, file_name, title):

    questions_columns = [
        ("Median Total Reads Mapped Correctly Across Samples", "Mapped_to_MT_Correct_Location"),
        ("Median Total Reads Mapped, but Not Correctly Across Samples", "Mapped_NOT_Correctly"),
        ("Median Reads Mapped Correctly (MQ > 30) Across Samples", "Mapped_to_MT_Correct_Location_MQ>30"),
        ("Median Total Reads Unmapped Across Samples", "Unmapped_Reads")
    ]

    fig, axes = plt.subplots(len(questions_columns), 1, figsize=(16, 25))
    fig.suptitle(title, fontsize=18, color='black')

    for ax, (label, column) in zip(axes, questions_columns):
        non_numeric_values = df[pd.to_numeric(df[column], errors='coerce').isna()][column]
        print(f"Non-numeric values in {column} column: {non_numeric_values.unique()}")
        assert all(pd.to_numeric(df[column], errors='coerce').notnull()), f"Non-numeric values found in {column} column."
        sns.barplot(
            x="Aligner_Name", 
            y=column, 
            hue="Damage_Type", 
            data=df, 
            ax=ax, 
            hue_order=['None', 'Mid', 'High', 'Single']
        )
        ax.set_title(label, fontsize=16, color='black')
        ax.set_xlabel("Alignment Algorithm", fontsize=14, color='black')
        ax.set_ylabel("Count of Reads", fontsize=14, color='black')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(file_name)

def create_latex_tables(df):
    metrics = [
        "Mapped_to_MT_Correct_Location",
        "Mapped_NOT_Correctly",
        "Mapped_to_MT_Correct_Location_MQ>30",
        "Unmapped_Reads"
    ]

    tables = {}

    # Abbreviations or short names for column headers
    abbreviations = {
        'SAFARI': 'SAF',
        'giraffe': 'GIR',
        'BWA-MEM': 'MEM',
        'BWA-aln (anc)': 'ALN-A',
        'BWA-aln': 'ALN',
        'BBMap': 'BBM',
        'SHRiMP': 'SHMP'
    }
    
    for metric in metrics:
        pivot_df = df.pivot_table(index='Damage_Type', columns='Aligner_Name', values=metric)
        pivot_df = pivot_df.rename(columns=abbreviations).round(2)  # Round to two decimals and rename columns
        latex_table = tabulate(pivot_df, tablefmt="latex_booktabs", headers="keys", showindex=True)
        caption = "\\caption{Placeholder caption}\n"
        label = f"\\label{{tab:{metric}}}\n"
        full_table = "\\begin{table}[ht]\n\\centering\n" + caption + label + latex_table + "\n\\end{table}"
        tables[metric] = full_table

    return tables


file_path = 'alignment_stats.csv'
df_new = pd.read_csv(file_path)
assert not df_new.empty, "The dataframe is empty."
processed_df = process_data(df_new)
create_new_plot(processed_df, "linear_benchmark.png", "Alignment Statistics Stratified by DNA Damage Type")

# Generate and save LaTeX tables
latex_tables = create_latex_tables(processed_df)

# Save the tables to individual files
for metric, table in latex_tables.items():
    with open(f"{metric}_table.tex", "w") as file:
        file.write(table)

# Print the tables to view them (optional)
for metric, table in latex_tables.items():
    print(table)

