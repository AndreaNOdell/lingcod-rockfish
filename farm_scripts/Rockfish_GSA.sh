#!/bin/bash -l

# setting name of job
#SBATCH --job-name=Rockfish_GSA

# setting home directory
#SBATCH -D /home/anodell/rockfish

# setting standard error output
#SBATCH -e /home/anodell/rockfish/slurm_log/sterror_%j.txt

# setting standard output
#SBATCH -o /home/anodell/rockfish/slurm_log/stdoutput_%j.txt

# setting medium priority
#SBATCH -p med

# setting the max time
#SBATCH -t 400:0:00

# mail alerts at beginning and end of job
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# send mail here
#SBATCH --mail-user=anodell@ucdavis.edu

# now we'll print out the contents of the R script to the standard output file
cat scripts/Rockfish_GSA.R
echo "ok now for the actual standard output"

# now running the actual script!

# load R
module load R

srun Rscript scripts/Rockfish_GSA.R
