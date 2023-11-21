#!/bin/bash
#SBATCH --job-name=hisat2-build      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=${USER}@ufl.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=20gb                    # Job memory request
#SBATCH --time=20:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date

module load edirect

cat $1 | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.gtf.gz/') ;
            wget "$url/$fname" ;
        done ;
    done

date
