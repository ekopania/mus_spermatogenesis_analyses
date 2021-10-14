#!/bin/bash
#PURPOSE: Run WGCNA
#
# Job name:
#SBATCH --job-name=WGCNA
#SBATCH --output=WGCNA-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
#SBATCH --mem=256G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate r4

#All data
#signed
Rscript 05_WGCNA.r all TRUE
#unsigned - NOT RECOMMENDED
#Rscript 05_WGCNA.r all FALSE

#Cell types separate
#signed
#Rscript 05_WGCNA.r LZ TRUE
#Rscript 05_WGCNA.r RS TRUE
#unsigned
#Rscript 05_WGCNA.r LZ FALSE
#Rscript 05_WGCNA.r RS FALSE

#Cell types separate, induced
#signed
#Rscript 05_WGCNA.r LZ TRUE
#Rscript 05_WGCNA.r RS TRUE
#unsigned
#Rscript 05_WGCNA.r LZ FALSE
#Rscript 05_WGCNA.r RS FALSE

echo "Done!"
