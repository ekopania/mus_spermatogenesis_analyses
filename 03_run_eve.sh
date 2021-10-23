#!/bin/bash
#PURPOSE: Run the EVE model (Rohlfs and Nielsen 2015)
#
# Job name:
#SBATCH --job-name=EVEmodel
#SBATCH --output=EVEmodel-%j.log
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

exp="/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/EVE_EXPRESSION_INPUTS/eve_expression_input_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm2.txt"
indiv="/mnt/beegfs/ek112884/mus_expression_analysis/eve_nindiv_input_four_species.txt"
tree="/mnt/beegfs/ek112884/mus_expression_analysis/eve_input_four_species.treefile"
ngene=$(head -1 ${exp})
run_name="LZexpressed_RPKM2"
subdir="RPKM2/"

echo "Running EVE command: /home/ek112884/software/EVE_release/EVEmodel -S -d ${exp} -i ${indiv} -t ${tree} -n ${ngene} -v 10 -f ${run_name} -p ${subdir}"

/home/ek112884/software/EVE_release/EVEmodel -S -d ${exp} -i ${indiv} -t ${tree} -n ${ngene} -v 10 -f ${run_name} -p ${subdir}

echo "Done!"

#FLAGS EXPLAINATION - copied from EVE README
#-S perform the EVE model expression divergence/diversity test on each gene using user provided expression data.
#-d the filename containing expression data
#-i the filename for the number of individuals per species
#-t the filename for the phylogeny file
#-n the number of genes in the dataset provided or simulated
#-v verbose describes how many updates are printed (0-100), values between 0 and 10 are recommended for standard use
#-f run specifier string is appended to result file names
#-p result path string is inserted between 'results/' and file names for result files
