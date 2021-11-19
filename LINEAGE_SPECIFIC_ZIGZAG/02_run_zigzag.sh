#!/bin/bash
#PURPOSE: Wrapper to run script 03_zigzag_plotLogExpDensity.r
#
# Job name:
#SBATCH --job-name=zigzag_plotLogExpDensity
#SBATCH --output=zigzag_plotLogExpDensity_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#ls *combined.tpm.txt | while read file; do
#ls erica_samples_*.combined.tpm.txt | while read file; do
#ls *spec.combined.tpm.txt | while read file; do
ls DOM_LZ.rmProbSample.spec.combined.tpm.txt | while read file; do
#	name=$(echo "${file}" | cut -d "." -f 1)
#	name=$(echo "${file}" | cut -d "." -f 1-2)
	name=$(echo "${file}" | cut -d "." -f 1-3)
	echo "${name}"
	Rscript 03_zigzag_plotLogExpDensity.r "${name}"
done
echo "Done!"
