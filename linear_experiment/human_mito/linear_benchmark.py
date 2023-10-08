import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def create_plot(df, file_name, title):
    # Check if all required columns are present
    required_columns = ['Damage_Type', 'Aligner_Name'] + [col for _, col in questions_columns]
    assert set(required_columns).issubset(df.columns), "Not all required columns are present in the dataframe."

    custom_order_df = df.copy()

    # Check if all values in Damage_Type are in the predefined order
    assert set(custom_order_df['Damage_Type']).issubset(damage_type_order), "Unexpected values found in Damage_Type column."
    
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].astype('category')
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].cat.reorder_categories(damage_type_order, ordered=True)
    custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].map(damage_type_rename)
    
    # Assert that all values in 'Aligner_Name' column are strings
    assert all(isinstance(name, str) for name in custom_order_df['Aligner_Name']), "Non-string values found in Aligner_Name column."
    
    custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].replace('safari', 'SAFARI')
    
    fig, axes = plt.subplots(len(journal_titles), 1, figsize=(16, 25))
    fig.suptitle(title, fontsize=18, color='black')
    
    for ax, (journal_title, column) in zip(axes, zip(journal_titles, [col for _, col in questions_columns])):
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
        ax.set_title(journal_title, fontsize=16, color='black')
        ax.set_xlabel("Alignment Algorithm", fontsize=14, color='black')
        ax.set_ylabel("Count of Reads", fontsize=14, color='black')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(file_name)

# Load the data
file_path = 'alignment_stats.csv'
df = pd.read_csv(file_path)

# Assert if the dataframe is empty
assert not df.empty, "The dataframe is empty."

# Define the questions and corresponding columns for plotting
questions_columns = [
    ("Number mapped to mt", "Mapped_to_MT"),
    ("Among those mapped, how many were mt and at the correct location+orientation?", "Mapped_to_MT_Correct_Location"),
    ("Among those mapped, how many were mt and at the correct location+orientation if MQ>30?", "Mapped_to_MT_Correct_Location_MQ>30"),
    ("Among those mapped, how many were not mt (i.e. bacteria+nuclear)?", "Mapped_NOT_to_MT"),
    ("Among those mapped, how many were not mt if MQ>30?", "Mapped_NOT_to_MT_MQ>30")
]

# Customize the Damage_Type levels and their order
damage_type_order = ['none', 'dmid', 'dhigh', 'single']
damage_type_rename = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single'}

# Define more descriptive journal-ready titles
journal_titles = [
    "Total Reads Mapped to mtDNA by Aligner",
    "Reads with Correct Location & Orientation in mtDNA",
    "Reads with Correct Location & Orientation in mtDNA (MQ > 30)",
    "Reads Mapped to Non-mtDNA Locations",
    "Reads Mapped to Non-mtDNA Locations (MQ > 30)"
]

# Create the original plot
create_plot(df, "linear_benchmark.png", "Alignment Statistics Stratified by DNA Damage Type")

