#!/bin/bash

# Directory containing the BAM files
BAM_DIR="alignments"

# Directory to save the output .prof files
OUTPUT_DIR="$BAM_DIR/profs"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Iterate through each BAM file in the BAM directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Get the base name of the BAM file (without path and extension)
    base_name=$(basename "$bam_file" .bam)

    # Run bam2prof on the BAM file and redirect the output to a .prof file in the output directory
    ./bam2prof/src/bam2prof "$bam_file" > "$OUTPUT_DIR/$base_name.prof"
done

