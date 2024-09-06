#!/bin/bash

#SBATCH --time=0-8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name="sc_Rscript"
#SBATCH --output="logs/slurm_%x_%j.out"
#SBATCH --error="logs/slurm_%x_%j.err"
#SBATCH --qos=short
#SBATCH --partition=c

mkdir -p logs/

ml load build-env/f2022
ml load r/4.1.2-foss-2021b
ml load r-bundle-bioconductor/3.14-foss-2021b-r-4.1.2

#Rscript script_regressOut.nCount_RNA.R
#Rscript script_regress_nCount_RNA.R
#Rscript script_run_bayesSpace.R
#Rscript script_snRNAseq_batchCorrection_animals.R
#Rscript script_scATAC_cisTopic_downsample.R
#Rscript script_scATAC_cisTopic.R
#Rscript run_LDA.R
Rscript run_RCTD.R
