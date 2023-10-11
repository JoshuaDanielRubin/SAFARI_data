import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def create_new_plot(df, file_name, title):
    # Check if all required columns are present
    required_columns = ['Damage_Type', 'Aligner_Name', 'Mapped_to_MT', 'Mapped_to_MT_Correct_Location', 'Mapped_to_MT_Correct_Location_MQ>30', 'Unmapped_Reads']
    assert set(required_columns).issubset(df.columns), "Not all required columns are present in the dataframe."

    custom_order_df = df.copy()

    # Create additional column for reads mapped but not correctly
    custom_order_df['Mapped_NOT_Correctly'] = custom_order_df['Mapped_to_MT'] - custom_order_df['Mapped_to_MT_Correct_Location']

    # Customize the Damage_Type levels and their order
    damage_type_order = ['none', 'dmid', 'dhigh', 'single']
    damage_type_rename = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single'}
    
    # Handle damage type order and renaming
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].astype('category')
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].cat.reorder_categories(damage_type_order, ordered=True)
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].map(damage_type_rename)

    # Handling Aligner Name
    aligner_rename = {
        'safari': 'SAFARI',
        'mem': 'BWA-MEM',
        'aln_anc': 'BWA-aln (anc)',
        'aln': 'BWA-aln',
        'bb': 'BBMap',
        'shrimp': 'SHRiMP'
    }
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].replace(aligner_rename)

    # Define a custom order for the Aligner Name and then sort the dataframe by this order
    custom_aligner_order = ['SAFARI', 'giraffe', 'BWA-MEM', 'BWA-aln (anc)', 'BWA-aln', 'BBMap', 'SHRiMP']
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].astype('category')
    custom_order_df['Aligner_Name'].cat.set_categories(custom_aligner_order, ordered=True, inplace=True)
    custom_order_df.sort_values('Aligner_Name', inplace=True)

    # Questions and corresponding columns for new plots
    questions_columns = [
        ("Median Total Reads Mapped Correctly Across Samples", "Mapped_to_MT_Correct_Location"),
        ("Median Total Reads Mapped, but Not Correctly Across Samples", "Mapped_NOT_Correctly"),
        ("Median Reads Mapped Correctly (MQ > 30) Across Samples", "Mapped_to_MT_Correct_Location_MQ>30"),
        ("Median Total Reads Unmapped Across Samples", "Unmapped_Reads")
    ]

    fig, axes = plt.subplots(len(questions_columns), 1, figsize=(16, 25))
    fig.suptitle(title, fontsize=18, color='black')

    for ax, (label, column) in zip(axes, questions_columns):
        # Check if values are numeric
        assert all(pd.to_numeric(custom_order_df[column], errors='coerce').notnull()), f"Non-numeric values found in {column} column."

        sns.barplot(
            x="Aligner_Name", 
            y=column, 
            hue="Damage_Type", 
            data=custom_order_df, 
            ax=ax, 
            hue_order=['None', 'Mid', 'High', 'Single']
        )
        ax.set_title(label, fontsize=16, color='black')
        ax.set_xlabel("Alignment Algorithm", fontsize=14, color='black')
        ax.set_ylabel("Count of Reads", fontsize=14, color='black')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(file_name)

    # Printing percent changes
    filtered_df = custom_order_df[custom_order_df['Aligner_Name'].isin(['giraffe', 'SAFARI'])]
    percent_changes = {}

    for _, column in questions_columns:
        pivot_df = filtered_df.pivot_table(index='Damage_Type', columns='Aligner_Name', values=column).reset_index()
        pivot_df['Percent Change'] = ((pivot_df['SAFARI'] - pivot_df['giraffe']) / pivot_df['giraffe']) * 100
        percent_changes[column] = pivot_df[['Damage_Type', 'Percent Change']]

    for metric, change_df in percent_changes.items():
        print(f"\nPercent change for {metric}:")
        for _, row in change_df.iterrows():
            damage_type, change = row['Damage_Type'], row['Percent Change']
            increase_or_decrease = "increase" if change > 0 else "decrease"
            print(f"{damage_type}: {abs(change):.2f}% {increase_or_decrease} from giraffe to SAFARI")


# Load the data
file_path = 'alignment_stats.csv'
df_new = pd.read_csv(file_path)

# Check if the dataframe is empty
assert not df_new.empty, "The dataframe is empty."

# Create the new plot
create_new_plot(df_new, "linear_benchmark.png", "Alignment Statistics Stratified by DNA Damage Type")

