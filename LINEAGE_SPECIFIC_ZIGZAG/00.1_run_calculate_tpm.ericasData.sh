#!/bin/bash
#PURPOSE: Wrapper to run 01.1_calculate_tpm.ericasData.r using slurm
#
# Job name:
#SBATCH --job-name=calculate_tpm_ericasData
#SBATCH --output=calculate_tpm_output.ericasData.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#Run R script to calculate tpm
Rscript 01.1_calculate_tpm.ericasData.r /mnt/beegfs/ek112884/RNAseqData_4cellType_fromErica/paired/CZIIhybrid_paired_suspenders_default_11Mar16_Q20_C_featureCounts.csv
echo "Done!"
