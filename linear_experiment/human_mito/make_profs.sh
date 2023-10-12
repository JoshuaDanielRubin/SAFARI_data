#!/bin/bash

# Directory containing the BAM files
BAM_DIR="alignments/with_md"

# Directory to save the output .prof files
OUTPUT_DIR="alignments/profs"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Iterate through each BAM file in the BAM directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Get the base name of the BAM file (without path and extension)
    base_name=$(basename "$bam_file" .bam)

    # Check if the base_name contains the word "single"
    if [[ "$base_name" == *"single"* ]]; then
        # Use the -single flag
        ./bam2prof/src/bam2prof -q -both -single -minl 20 "$bam_file" > "$OUTPUT_DIR/$base_name.prof"
    else
        # Use the -double flag
        ./bam2prof/src/bam2prof -q -paired -both -double -minl 20 "$bam_file" > "$OUTPUT_DIR/$base_name.prof"
    fi
done

