#!/bin/bash

# Define number of cores to use
num_cores=20

# Create full directory if it doesn't already exist
mkdir -p full

# Use xargs to run wget in parallel across multiple cores
cat ENA_ids.txt | nice xargs -n 1 -P $num_cores wget -P full

