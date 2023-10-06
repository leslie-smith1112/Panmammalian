#!/bin/bash
#SBATCH --job-name=merge_human_cancers      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=300gb                    # Job memory request
#SBATCH --time=080:00:00              # Time limit hrs:min:sec
#SBATCH --output=merge_human_cancers%j.log     # Standard output and error log
#SBATCH --error=merge_human_cancers%j.log #error log

pwd; hostname; date

ml R

Rscript --max-ppsize=500000 /home/leslie.smith1/blue_kgraim/leslie.smith1/pan/merge_pancancer_data.R

date
