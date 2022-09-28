#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name="sc_pyscript"
#SBATCH --output="logs/slurm-%x_%j.out"
#SBATCH --error="logs/slurm-%x_%j.err"
#SBATCH --qos=medium
#SBATCH --partition=c

mkdir -p logs/

python script_cell2location_subtypes.py

