#!/bin/bash
#PURPOSE: Run WGCNA data checks
#
# Job name:
#SBATCH --job-name=WGCNA_dataChecks
#SBATCH --output=WGCNA_dataChecks-%j.log
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
#Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/03_WGCNA_dataChecks.r all

#Cell types separate
#Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/03_WGCNA_dataChecks.r LZ
Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/03_WGCNA_dataChecks.r RS

#Cell types separate, induced only
#Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/03_WGCNA_dataChecks.r LZind
#Rscript /mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/mus_spermatogenesis_analyses/WGCNA_scripts/03_WGCNA_dataChecks.r RSind

echo "Done!"
