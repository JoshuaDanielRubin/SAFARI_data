#!/bin/bash

# Input parameter
bam_dir=$1  # Path to directory containing BAM files

# Create directories for the subsampled BAM files at different depths
#mkdir -p subsampled/0.25x
#mkdir -p subsampled/1x
#mkdir -p subsampled/2x
#mkdir -p subsampled/0.5x

# Function to subsample a BAM file at a given depth
subsample_bam() {
  bam_file=$1
  chr=$2
  target_depth=$3
  replicate=$4

  # Calculate current depth of coverage on the mitochondria
  depth=$(samtools depth -r $chr $bam_file | awk '{sum+=$3} END {print sum/NR}')

  # Generate a random seed
  seed=$(shuf -i 1-100000 -n 1)

  # Calculate the fraction of reads to subsample using awk for better precision
  fraction=$(awk -v target=$target_depth -v depth=$depth 'BEGIN {printf "%.6f", target/depth}')

  # Construct the argument for samtools view -bs in the correct SEED.FRACTION format
  sampling_argument="${seed}${fraction#0}"

  # Subsample the BAM file
  filename=$(basename $bam_file)
  out_dir="subsampled_reps/${target_depth}x"
  out_file="${filename%%.bam}_replicate_${replicate}.bam"
  echo "Seed: $seed, Fraction: $fraction, Sampling Argument: $sampling_argument"
  nice -19 samtools view -bs $sampling_argument $bam_file -@ 60 > $out_dir/$out_file

  # Calculate the depth of coverage of the subsampled BAM file
  subsampled_depth=$(samtools depth $out_dir/$out_file | awk '{sum+=$3} END {if (NR > 0) print sum/NR; else print "0"}')

  echo "Depth of coverage of $bam_file on the mitochondria is $depth"
  echo "Subsampled BAM file to a depth of $subsampled_depth"
  echo "Subsampled BAM file saved as $out_dir/$out_file"
}

# Export the subsample_bam function so it can be used by GNU Parallel
export -f subsample_bam

# Loop through all BAM files in the input directory
for bam_file in ${bam_dir}/*.bam; do

  # Check if mitochondrial chromosome is called "chrM" or "MT"
  if samtools idxstats $bam_file | cut -f1 | grep -q "chrM"; then
    chr="chrM"
  elif samtools idxstats $bam_file | cut -f1 | grep -q "MT"; then
    chr="MT"
  else
    echo "Mitochondrial chromosome not found in $bam_file"
    continue
  fi

  # Subsample the BAM file using GNU Parallel
  filename=$(basename $bam_file)
  echo "Subsampling $filename..."

  # Subsample at target depth of 0.25x, for 5 replicates
  #for replicate in {1..5}; do
  #  nice -19 parallel -j 20 subsample_bam {} $chr 0.25 $replicate ::: $bam_file
  #done

  # Subsample at target depth of 1x, for 5 replicates
  #for replicate in {1..5}; do
  #  nice -19 parallel -j 60 subsample_bam {} $chr 1 $replicate ::: $bam_file
  #done

  # Subsample at target depth of 2x
   #for replicate in {1..5}; do
   #  nice -19 parallel -j 60 subsample_bam {} $chr 2 $replicate ::: $bam_file
   #done

  # Subsample at target depth of 0.5x
   #for replicate in {1..5}; do
   #  nice -19 parallel -j 55 subsample_bam {} $chr 0.5 $replicate ::: $bam_file
   #done

done

