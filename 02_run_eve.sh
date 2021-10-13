#!/bin/bash
#PURPOSE: Run the EVE model (Rohlfs and Nielsen 2015)
#
# Job name:
#SBATCH --job-name=EVEmodel
#SBATCH --output=EVEmodel-%j.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-3 # run on a specific node
#
## Command(s) to run:

/home/ek112884/software/EVE_release/EVEmodel <flags>
