import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load the data
file_path = 'alignment_stats.csv'
df = pd.read_csv(file_path)

# Define the questions and corresponding columns for plotting
questions_columns = [
    ("Number mapped to mt", "Mapped_to_MT"),
    ("Among those mapped, how many were mt and at the correct location+orientation?", "Mapped_to_MT_Correct_Location"),
    ("Among those mapped, how many were mt and at the correct location+orientation if MQ>30?", "Mapped_to_MT_Correct_Location_MQ>30"),
    ("Among those mapped, how many were not mt (i.e. bacteria+nuclear)?", "Mapped_NOT_to_MT"),
    ("Among those mapped, how many were not mt if MQ>30?", "Mapped_NOT_to_MT_MQ>30")
]

# Customize the Damage_Type levels and their order
damage_type_order = ['dnone', 'ddmid', 'ddhigh', 'dsingle']
custom_order_df = df.copy()
custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].astype('category')
custom_order_df['Damage_Type'].cat.reorder_categories(damage_type_order, ordered=True, inplace=True)

# Rename the Damage_Type levels for better readability
damage_type_rename = {'dnone': 'None', 'ddmid': 'Mid', 'ddhigh': 'High', 'dsingle': 'Single'}
custom_order_df['Damage_Type'] = custom_order_df['Damage_Type'].map(damage_type_rename)

# Update Aligner_Name to make 'safari' uppercase ('SAFARI')
custom_order_df['Aligner_Name'] = custom_order_df['Aligner_Name'].replace('safari', 'SAFARI')

# Define more descriptive journal-ready titles
journal_titles = [
    "Total Reads Mapped to mtDNA by Aligner",
    "Reads with Correct Location & Orientation in mtDNA",
    "Reads with Correct Location & Orientation in mtDNA (MQ > 30)",
    "Reads Mapped to Non-mtDNA Locations",
    "Reads Mapped to Non-mtDNA Locations (MQ > 30)"
]

# Initialize the figure
fig, axes = plt.subplots(len(journal_titles), 1, figsize=(16, 25))
fig.suptitle("Alignment Statistics Stratified by DNA Damage Type", fontsize=18, color='black')

# Loop over each question and plot the data
for ax, (journal_title, column) in zip(axes, zip(journal_titles, [col for _, col in questions_columns])):
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

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save the plot to a file
plt.savefig("linear_benchmark.png")

