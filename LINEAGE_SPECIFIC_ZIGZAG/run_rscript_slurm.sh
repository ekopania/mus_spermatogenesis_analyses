#!/bin/bash
#PURPOSE: Generic slurm wrapper for running an r script. At the bottom of this script under #Commands to run, type:
#Rscript <name_of_script.r> <command line args>
#
# Job name:
#SBATCH --job-name=rscript
#SBATCH --output=rscript_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#Rscript calculate_rpkm.r
#Rscript calculate_rpkm.ericasdata.r
#Rscript analyze_zigzag_output.WDIS.r
Rscript 08_analyze_zigzag_output.spec.r
