#!/bin/bash

cat <(zcat $1) <(/home/ctools/interleave_fastq/interleave_fastq.sh $2  $3)

