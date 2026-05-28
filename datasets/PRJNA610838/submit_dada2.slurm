#!/bin/bash
#SBATCH --job-name=PRJNA610838_dada2
#SBATCH --account=gts-lspeare3
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH -t 00:30:00
#SBATCH -p cpu-small
#SBATCH -o PRJNA610838_dada2_%j.log
#SBATCH -e PRJNA610838_dada2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rthomas376@gatech.edu

# Load R module
module load r/4.3.3

# Run your R script
Rscript PRJNA610838_dada2_batch_job.R
