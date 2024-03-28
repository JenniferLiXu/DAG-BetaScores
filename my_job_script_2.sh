#!/bin/bash
#SBATCH --job-name=my_job       # Name of the job
#SBATCH --output=my_job_%j.out  # Standard output and error log (%j expands to jobID)
#SBATCH --error=my_job_%j.err   # Standard error log
#SBATCH --time=05:00:00         # Time limit hrs:min:sec
#SBATCH --partition=main        # Partition where the job should be submitted
#SBATCH --ntasks=1              # Run on a single CPU
#SBATCH --mem-per-cpu=4G                # Memory limit
module load R
Rscript test_singleloop_2.R


