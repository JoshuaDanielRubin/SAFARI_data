#!/usr/bin/python
import pysam
import sys

# Function to check if two intervals intersect
def intersects(left1, right1, left2, right2):
    return not (right1 < left2 or left1 > right2)

def main(bam_file_path):
    # Suppress SAM-format errors
    pysam.set_verbosity(0)

    # Open the BAM file for reading
    bamInputFile = pysam.Samfile(bam_file_path, "rb")

    # Initialize counters
    total = 0
    mapped_to_mt = 0
    mapped_to_mt_correct_location = 0
    mapped_to_mt_correct_location_MQ_gt_30 = 0
    not_mapped_to_mt = 0
    not_mapped_to_mt_MQ_gt_30 = 0

    # Loop through each read in the BAM file
    for read in bamInputFile:
        total += 1  # Increment total read counter

        # Only consider mapped reads
        if not read.is_unmapped:
            # Check if read is from mitochondrial genome based on query name
            if read.query_name.startswith("generation"):
                mapped_to_mt += 1
                
                # Check for correct location and orientation
                fields = read.query_name.split(":")
                intcs = intersects(
                    int(fields[2]), int(fields[3]),
                    read.reference_start - 50, read.reference_start + read.reference_length + 50
                )
                if intcs:
                    mapped_to_mt_correct_location += 1
                    
                    # Check for MQ > 30
                    if read.mapping_quality > 30:
                        mapped_to_mt_correct_location_MQ_gt_30 += 1
            else:
                not_mapped_to_mt += 1
                
                # Check for MQ > 30
                if read.mapping_quality > 30:
                    not_mapped_to_mt_MQ_gt_30 += 1

    # Close the BAM file
    bamInputFile.close()

    # Assert that all reads are accounted for
    #assert total == mapped_to_mt + not_mapped_to_mt, "Not all reads are accounted for!"

    # Output statistics in table form
    print(f"Total reads: {total}")
    print(f"-----------------------------------------")
    print(f"Mapped to MT: {mapped_to_mt}")
    print(f"Mapped to MT (Correct Location): {mapped_to_mt_correct_location}")
    print(f"Mapped to MT (Correct Location, MQ>30): {mapped_to_mt_correct_location_MQ_gt_30}")
    print(f"NOT Mapped to MT: {not_mapped_to_mt}")
    print(f"NOT Mapped to MT (MQ>30): {not_mapped_to_mt_MQ_gt_30}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <bam_file_path>")
    else:
        main(sys.argv[1])

