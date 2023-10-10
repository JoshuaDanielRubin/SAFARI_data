# SAFARI Experiment Data

This repository contains the scripts and data necessary to reproduce our benchmarking experiment results. Please note that you may need to modify some paths to reproduce the results on your machine.

## Scripts

- **Euka Experiment**: 
  - Script: `./run_euka.sh`
  - Description: Runs the euka experiment.

- **HaploCart Experiment**: 
  - Script: `./run_hc.sh`
  - Description: Runs the HaploCart experiment.

- **Index Creation**: 
  - Script: `./index.sh`
  - Description: Creates the minimizer and RYmer index files on disk.

## linear_experiment/human_mito Directory

- **Evaluate**: 
  - Script: `evaluate.py`
  - Description: Computes summary statistics from the `alignment_stats.csv` file.

- **Compare Damage Estimate**: 
  - Script: `compare_damage_est.py`
  - Description: Performs a cross-comparison of damage rate estimates from bam2prof.

- **Interleave FASTQ**: 
  - Script: `interleave_fastq.sh`
  - Description: Interleaves FASTQ files.

- **Populate MD Field**: 
  - Script: `make_MD.sh`
  - Description: Populates the MD field in the BAM files based on the reference, as required by bam2prof.

- **Create Damage Profiles**: 
  - Script: `make_profs.sh`
  - Description: Runs bam2prof to create the damage profiles (i.e., estimated nucleotide misincorporation matrices).

- **Generate simulations for linear aligner benchmarking experiment**:
   - Script: `Snakefile`
   - Description: Generates simulations

- **Compute statistics from BAM file**
   - Script: `parseBamMito.py`
   - Description: Compute statistics from BAM file

## spurious_model Directory

- **Estimate parameters for spurious alignment model**
   - Script: `main.py`
   - Description: Estimates parameters of the power law distribution for our spurious alignment model using a RYmer and minimizer index.
