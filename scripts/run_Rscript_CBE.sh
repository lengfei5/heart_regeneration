#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --job-name="sc_Rscript"
#SBATCH --output="logs/slurm-%x_%j.out"
#SBATCH --error="logs/slurm-%x_%j.err"
#SBATCH --qos=medium
#SBATCH --partition=c

mkdir -p logs/

ml load build-env/f2022
ml load r/4.1.2-foss-2021b

#Rscript script_regressOut.nCount_RNA.R
Rscript script_regress_nCount_RNA.R
