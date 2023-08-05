#!/bin/bash

# Ensure two arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: ./kmer_stats.sh <fasta_file> <kmer_size>"
    exit 1
fi

# Get the arguments
fasta_file=$1
kmer_size=$2

# Define the output file names
jf_file="output.jf"
txt_file="kmers_${kmer_size}.txt"
gz_file="${txt_file}.gz"

# Run jellyfish count
jellyfish count -m ${kmer_size} -s 100M -t 60 -o ${jf_file} ${fasta_file}

# Run jellyfish dump
jellyfish dump -c -t ${jf_file} > ${txt_file}

# Gzip the text file
gzip ${txt_file}

# Remove the .jf file
rm ${jf_file}

# Print completion message
echo "K-mer counting complete. Results are in ${gz_file}."

