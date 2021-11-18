#!/bin/bash
#PURPOSE: Merge and log transform normalized expression data for input into WGCNA; make table of "traits" (cross type and cell type)
#
# Job name:
#SBATCH --job-name=WGCNA_preprocess
#SBATCH --output=WGCNA_preprocess-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=4G #Not sure if I should mess with these...
#SBATCH --mem=16G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/01_process_exp_data.r

echo "Done!"
