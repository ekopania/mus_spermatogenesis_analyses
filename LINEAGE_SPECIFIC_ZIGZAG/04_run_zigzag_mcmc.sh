#!/bin/bash
#PURPOSE: Wrapper to run script 05_zigzag_mcmc.r
#
# Job name:
#SBATCH --job-name=zigzag_mcmc
#SBATCH --output=zigzag_mcmc_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=96000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#ls *combined.tpm.txt | while read file; do
ls PAHNew_RS.combined.tpm.txt | while read file; do
#ls *_RS.spec.combined.tpm.txt | while read file; do
	name=$(echo "${file}" | cut -d "." -f 1)
#	name=$(echo "${file}" | cut -d "." -f 1-2)
#	name=$(echo "${file}" | cut -d "." -f 1-3)
	echo "${name}"
	Rscript 05_zigzag_mcmc.r "${name}"
done
echo "Done!"
