#!/bin/bash
#SBATCH --job-name=parallel_job      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rveernapu@ufl.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=32            # Number of CPU cores per task
#SBATCH --mem=128gb                    # Job memory request
#SBATCH --time=25:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date

cd /blue/kgraim/common/genome

wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.hal

date

