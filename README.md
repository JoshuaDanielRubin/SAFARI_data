# SAFARI
Sensitive Alignments from a RYmer Index

This repository contains a modified version of vg (https://github.com/vgteam/vg) in which the subcommand vg giraffe has been modified to be more suited for ancient DNA samples.

The repository also contains the scripts and data necessary to reproduce our benchmarking experiment results. Note that some paths may need to be modified to reproduce the results on your machine.

The scripts are as follows:

   - ./run_euka.sh runs the euka experiment.
   - ./run_giraffe.sh runs the giraffe experiment (i.e. direct comparison of corrected vs. uncorrected in a pangenome context)
   - ./run_hc.sh runs the HaploCart experiment
   - ./index.sh creates the minimizer and RYmer index files on disk.
