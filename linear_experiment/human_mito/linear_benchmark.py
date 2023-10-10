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

    # Handle damage type order and renaming
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].astype('category')
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].cat.reorder_categories(damage_type_order, ordered=True)
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].map(damage_type_rename)

    # Handling Aligner Name
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].replace('safari', 'SAFARI')

    # Calculate the mean total number of reads
    mean_total_reads = df[['Mapped_to_MT', 'Mapped_to_MT_Correct_Location', 'Mapped_to_MT_Correct_Location_MQ>30', 'Unmapped_Reads']].sum(axis=1).mean()

    # Define the questions and corresponding columns for new plots
    questions_columns = [
        ("Total Reads Mapped Correctly", "Mapped_to_MT_Correct_Location"),
        ("Total Reads Mapped, but Not Correctly", "Mapped_NOT_Correctly"),
        ("Reads Mapped Correctly (MQ > 30)", "Mapped_to_MT_Correct_Location_MQ>30"),
        ("Total Reads Unmapped", "Unmapped_Reads")
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
        ax.set_ylim(0, mean_total_reads)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(file_name)

# Load the new data
file_path = 'alignment_stats.csv'
df_new = pd.read_csv(file_path)

# Check if the dataframe is empty
assert not df_new.empty, "The dataframe is empty."

# Customize the Damage_Type levels and their order
damage_type_order = ['none', 'dmid', 'dhigh', 'single']
damage_type_rename = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single'}

# Create the new plot
create_new_plot(df_new, "linear_benchmark.png", "Alignment Statistics Stratified by DNA Damage Type")

