#!/bin/sh
#SBATCH --job-name=DEXSeqCountMatrix   # Job name
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu # Where to send mail	
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=5:00:00             # Time limit hrs:min:sec
# This is an example script that combines array tasks with
# bash loops to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users. 
pwd; hostname; date

module load R
#arguments: dataset, mammal

DATASET=$1
MAMMAL=canis_familiaris

echo "dataset being processed is $DATASET, mammal is $MAMMAL"

Rscript DEXSeq.R -d ${DATASET} -m ${MAMMAL}


date
