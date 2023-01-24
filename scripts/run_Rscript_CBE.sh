#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --job-name="sc_Rscript"
#SBATCH --output="logs/slurm-%x_%j.out"
#SBATCH --error="logs/slurm-%x_%j.err"
#SBATCH --qos=medium
#SBATCH --partition=m

mkdir -p logs/

ml load build-env/f2022
ml load r/4.1.2-foss-2021b
ml load r-bundle-bioconductor/3.14-foss-2021b-r-4.1.2

#Rscript script_regressOut.nCount_RNA.R
#Rscript script_regress_nCount_RNA.R
#Rscript script_run_bayesSpace.R
#Rscript script_snRNAseq_batchCorrection_animals.R
#Rscript script_scATAC_cisTopic_downsample.R
Rscript script_scATAC_cisTopic.R
