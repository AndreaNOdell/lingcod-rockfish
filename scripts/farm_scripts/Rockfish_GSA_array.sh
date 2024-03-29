#!/bin/bash -l

# setting name of job
#SBATCH --job-name=Rockfish_GSA_array

# setting home directory
#SBATCH -D /home/anodell/rockfish

# Setting an array
#SBATCH --array=0-19

# setting standard error output
#SBATCH -e /home/anodell/rockfish/slurm_log/sterror_Rockfish_GSA_array_%A_%a.txt

# setting standard output
#SBATCH -o /home/anodell/rockfish/slurm_log/stdoutput_Rockfish_GSA_array_%A_%a.txt

#SBATCH --account=baskettgrp

# setting medium priority
#SBATCH -p med2

# setting the max time
#SBATCH -t 700:0:00

# mail alerts at beginning and end of job
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# send mail here
#SBATCH --mail-user=anodell@ucdavis.edu

# now we'll print out the contents of the R script to the standard output file
cat scripts/Rockfish_GSA_array.R
echo "ok now for the actual standard output"

# now running the actual script!

# load R
module load R

srun Rscript scripts/Rockfish_GSA_array.R ${SLURM_ARRAY_TASK_ID}
