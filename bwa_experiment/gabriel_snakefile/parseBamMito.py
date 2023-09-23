#!/usr/bin/python

import pysam
import sys

pysam.set_verbosity(0)
bamInputFile = pysam.Samfile(sys.argv[1], "rb")

mapped = 0
unmapped = 0
correctlocationmapped = 0
numt = 0
total = 0

def intersects(left1, right1, left2, right2):
    return not (right1 < left2 or left1 > right2)

for read in bamInputFile:
    total += 1

    if read.is_unmapped:
        unmapped += 1
    else:
        mapped += 1  # Now counting all mapped reads, not just those with a specific query name prefix
        
        if read.query_name[0:10] == "generation":
            fields = read.query_name.split(":")
            intcs = intersects(
                int(fields[2]), int(fields[3]),
                read.reference_start-50, read.reference_start+read.reference_length+50
            )
            if intcs:
                correctlocationmapped += 1
        else:
            numt += 1

bamInputFile.close()

print("Total:\t" + str(total))
print("unmapped:\t" + str(unmapped))
print("unmapped%:\t" + str(100.0 * (unmapped / total) if total != 0 else 0.0))
print("mapped:\t" + str(mapped))
print("mapped%:\t" + str(100.0 * (mapped / total) if total != 0 else 0.0))
print("correctmap:\t" + str(correctlocationmapped))
print("correctmap%:\t" + str(100.0 * (correctlocationmapped / mapped) if mapped != 0 else 0.0))
print("numt:\t" + str(numt))
print("numt%:\t" + str(100.0 * (numt / total) if total != 0 else 0.0))

