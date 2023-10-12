import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def create_new_plot(df, file_name, title):
    required_columns = ['Damage_Type', 'Aligner_Name', 'Mapped_to_MT', 'Mapped_to_MT_Correct_Location', 'Mapped_to_MT_Correct_Location_MQ>30', 'Unmapped_Reads']
    assert set(required_columns).issubset(df.columns), "Required columns missing."
    
    df['Mapped_NOT_Correctly'] = df['Mapped_to_MT'] - df['Mapped_to_MT_Correct_Location']
    df['Damage_Type'] = df['Damage_Type'].astype('category').cat.reorder_categories(damage_type_order, ordered=True).map(damage_type_rename)
    df['Aligner_Name'] = df['Aligner_Name'].replace('safari', 'SAFARI')
    questions_columns = [("Median Total Reads Mapped Correctly", "Mapped_to_MT_Correct_Location"), ("Median Total Reads Mapped, but Not Correctly", "Mapped_NOT_Correctly"), ("Median Reads Mapped Correctly (MQ > 30)", "Mapped_to_MT_Correct_Location_MQ>30"), ("Median Total Reads Unmapped", "Unmapped_Reads")]

    fig, axes = plt.subplots(len(questions_columns), 1, figsize=(16, 25))
    fig.suptitle(title, fontsize=18, color='black')

    for ax, (label, column) in zip(axes, questions_columns):
        sns.barplot(x="Aligner_Name", y=column, hue="Damage_Type", data=df, ax=ax, hue_order=['None', 'Mid', 'High', 'Single'])
        ax.set_title(label, fontsize=16)
        ax.set_xlabel("Alignment Algorithm", fontsize=14)
        ax.set_ylabel("Median Count of Reads", fontsize=14)
        
        # Compute Median values for each aligner by damage type
        table_df = df.groupby(['Aligner_Name', 'Damage_Type'])[column].median().unstack()
        
        # Print Latex table
        table = table_df.to_latex()
        print(f"\nLaTeX Table for {label}:\n", table)

        # Calculate and print percent increase or decrease from "giraffe" to "SAFARI" for each damage type
        for damage in damage_type_rename.values():
            giraffe_val = table_df.loc["giraffe", damage]
            safari_val = table_df.loc["SAFARI", damage]
            if giraffe_val:
                percent_change = ((safari_val - giraffe_val) / giraffe_val) * 100
                print(f"Percent change from 'giraffe' to 'SAFARI' for {label} with damage type '{damage}': {percent_change:.2f}%")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(file_name)

# Load data and check
file_path = 'alignment_stats.csv'
df_new = pd.read_csv(file_path)
assert not df_new.empty, "Dataframe empty."

# Customize the Damage_Type levels and order
damage_type_order = ['none', 'dmid', 'dhigh', 'single']
damage_type_rename = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single'}
create_new_plot(df_new, "linear_benchmark.png", "Alignment Statistics Stratified by DNA Damage Type")

