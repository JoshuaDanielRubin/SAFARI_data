import pandas as pd
import gzip
from Levenshtein import distance as levenshtein_distance

# Function to read a zipped FASTA file and return the consensus sequence
def read_fasta(file_path):
    with gzip.open(file_path, 'rt') as f:
        # Read the sequence lines and concatenate them
        lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return sequence

# Function to calculate mean edit distances for uncorrected and corrected sequences
def calculate_mean_edit_distance(table_path, fasta_dir):
    # Read the table
    df = pd.read_csv(table_path, sep=' & ', engine='python', header=None)
    
    # Extract relevant columns
    ground_truth_haplogroups = df[2].tolist()
    corrected_haplogroups = df[4].tolist()
    uncorrected_haplogroups = df[5].tolist()

    edit_distances_uncorrected = []
    edit_distances_corrected = []
    for ground_truth, corrected, uncorrected in zip(ground_truth_haplogroups, corrected_haplogroups, uncorrected_haplogroups):
        # Construct file paths and sanitize them
        ground_truth_path = f"{fasta_dir}/{ground_truth}.fasta.gz".replace(" ", "").replace("\t", "")
        corrected_path = f"{fasta_dir}/{corrected}.fasta.gz".replace(" ", "").replace("\t", "")
        uncorrected_path = f"{fasta_dir}/{uncorrected}.fasta.gz".replace(" ", "").replace("\t", "")

        # Read the consensus sequences
        ground_truth_seq = read_fasta(ground_truth_path)
        corrected_seq = read_fasta(corrected_path)
        uncorrected_seq = read_fasta(uncorrected_path)

        # Calculate and store the edit distances
        distance_uncorrected = levenshtein_distance(ground_truth_seq, uncorrected_seq)
        distance_corrected = levenshtein_distance(ground_truth_seq, corrected_seq)

        edit_distances_uncorrected.append(distance_uncorrected)
        edit_distances_corrected.append(distance_corrected)

    # Calculate means
    mean_distance_uncorrected = pd.Series(edit_distances_uncorrected).mean()
    mean_distance_corrected = pd.Series(edit_distances_corrected).mean()

    return mean_distance_uncorrected, mean_distance_corrected

# Example usage
table_path = 'raw_table.txt'
fasta_dir = '/home/projects/mito_haplotype/vgan/data/synthetic_fastas'
mean_distance_uncorrected, mean_distance_corrected = calculate_mean_edit_distance(table_path, fasta_dir)
print(f"mean Edit Distance (Uncorrected vs Ground Truth): {mean_distance_uncorrected}")
print(f"mean Edit Distance (Corrected vs Ground Truth): {mean_distance_corrected}")

