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
    unmapped_reads = 0  # Initialize counter for unmapped reads

    # Loop through each read in the BAM file
    for read in bamInputFile:
        total += 1  # Increment total read counter

        # Check for unmapped reads
        if read.is_unmapped:
            unmapped_reads += 1
        else:
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

    # Assertions
    assert total == (mapped_to_mt + not_mapped_to_mt + unmapped_reads), "Mismatch in total mapped and unmapped reads!"
    assert mapped_to_mt == (mapped_to_mt_correct_location + (mapped_to_mt - mapped_to_mt_correct_location)), "Mismatch in MT mapped reads!"
    assert mapped_to_mt_correct_location_MQ_gt_30 <= mapped_to_mt_correct_location, "More reads with MQ>30 than reads in the correct location!"
    assert not_mapped_to_mt_MQ_gt_30 <= not_mapped_to_mt, "More non-MT reads with MQ>30 than total non-MT reads!"
    assert (mapped_to_mt + not_mapped_to_mt) == (total - unmapped_reads), "Mismatch in MT vs non-MT reads!"

    # Output statistics in table form
    print(f"Total reads: {total}")
    print(f"-----------------------------------------")
    print(f"Unmapped reads: {unmapped_reads}")
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

