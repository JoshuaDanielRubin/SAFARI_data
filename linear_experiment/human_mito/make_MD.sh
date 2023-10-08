#!/bin/bash

# Define the directories and reference FASTA files
source_dir="alignments"
target_dir="alignments/with_md"
reference_fasta1="simulations/gen_0.fa"
reference_fasta2="rCRS.fa"

# Ensure samtools is available
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found. Please install samtools."
    exit
fi

# Create the target directory if it doesn't exist
mkdir -p "$target_dir"

# Iterate through all BAM files in the source directory
for bam_file in "$source_dir"/*.bam; do
    # Extract the filename without the directory
    filename=$(basename "$bam_file")

    # Construct the new file path in the target directory
    new_file="$target_dir/$filename"

    # Determine the correct reference file by inspecting the header of the BAM file
    if samtools view -H "$bam_file" | grep -q "SN:generation_0"; then
        reference_fasta="$reference_fasta1"
    else
        reference_fasta="$reference_fasta2"
    fi

    # Recalculate the MD tags and save the new BAM file
    if ! samtools calmd -b "$bam_file" "$reference_fasta" > "$new_file"
    then
        # Output an error message if the command fails
        echo "Error processing $bam_file with reference file."
        continue  # Skip to the next iteration
    fi

    # Output progress
    echo "Processed $bam_file --> $new_file"
done

echo "Processing complete."

