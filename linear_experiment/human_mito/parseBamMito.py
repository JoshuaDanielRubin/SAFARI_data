#!/usr/bin/python
import pysam
import sys

def main(bam_file_path):
    # Suppress SAM-format errors
    pysam.set_verbosity(0)

    # Open the BAM file for reading
    bamInputFile = pysam.Samfile(bam_file_path, "rb")

    # Initialize counters for TP, FP, TN, FN and their counterparts with MQ > 30
    TP, FP, TN, FN = 0, 0, 0, 0
    TP_MQ30, FP_MQ30, TN_MQ30, FN_MQ30 = 0, 0, 0, 0

    # Loop through each read in the BAM file
    for read in bamInputFile:
        # Determine if read came from MT
        is_from_mt = read.query_name.startswith("generation")

        # Check if read is mapped to MT
        is_mapped_to_mt = not read.is_unmapped and (read.reference_name == "generation_0" or read.reference_name == "H2a2a1")

        # Compute TP, FP, TN, and FN
        if is_from_mt and is_mapped_to_mt:
            TP += 1
            if read.mapping_quality > 30:
                TP_MQ30 += 1
        elif not is_from_mt and is_mapped_to_mt:
            FP += 1
            if read.mapping_quality > 30:
                FP_MQ30 += 1
        elif is_from_mt and not is_mapped_to_mt:
            FN += 1
            FN_MQ30 += 1  # Unmapped reads from MT 
        elif not is_from_mt and not is_mapped_to_mt:
            TN += 1
            TN_MQ30 += 1  # Unmapped reads not from MT

    # Close the BAM file
    bamInputFile.close()

    # Output the computed values
    print(f"True Positives (TP): {TP}")
    print(f"False Positives (FP): {FP}")
    print(f"True Negatives (TN): {TN}")
    print(f"False Negatives (FN): {FN}")
    print("\nFor reads with MQ > 30:")
    print(f"True Positives (TP_MQ30): {TP_MQ30}")
    print(f"False Positives (FP_MQ30): {FP_MQ30}")
    print(f"True Negatives (TN_MQ30): {TN_MQ30}")
    print(f"False Negatives (FN_MQ30): {FN_MQ30}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <bam_file_path>")
    else:
        main(sys.argv[1])

