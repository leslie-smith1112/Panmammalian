#!/bin/bash
#SBATCH --job-name=parallel_job      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=${USER}@ufl.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=16gb                    # Job memory request
#SBATCH --time=01:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date

echo ${USER}

# Runs very fast - not much memory is needed
cd /blue/kgraim/common/genome

module load samtools
samtools faidx $1

module load picard

picard CreateSequenceDictionary -R $1  
#      O=reference.dict

date

