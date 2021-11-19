#!/bin/bash
#PURPOSE: extract omegas for gene lists based on ensembl pairwise 1:1 ortho dataset
#
# Job name:
#SBATCH --job-name=extract_omegas
#SBATCH --output=extract_omegas-%j.log
##SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=8G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro

myDir="MULTI_MAP"
rpkm=10
induced_cutoff=10

echo "Extracting omegas for ${myDir} at an induced cutoff of ${induced_cutoff}"

rm LZ_expressed_no_paml.ensemblOrthos.${myDir}.rpkm${rpkm}.txt
rm omega_list_LZ.ensemblOrthos.${myDir}.rpkm${rpkm}.txt
cat /mnt/beegfs/ek112884/mus_expression_analysis/${myDir}/gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm${rpkm}.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_expressed_no_paml.ensemblOrthos.${myDir}.rpkm${rpkm}.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZ.ensemblOrthos.${myDir}.rpkm${rpkm}.txt; done

rm RS_expressed_no_paml.ensemblOrthos.${myDir}.rpkm${rpkm}.txt
rm omega_list_RS.ensemblOrthos.${myDir}.rpkm${rpkm}.txt
cat /mnt/beegfs/ek112884/mus_expression_analysis/${myDir}/gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.rpkm${rpkm}.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_expressed_no_paml.ensemblOrthos.${myDir}.rpkm${rpkm}.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RS.ensemblOrthos.${myDir}.rpkm${rpkm}.txt; done

#rm LZ_induced${induced_cutoff}_no_paml.ensemblOrthos.${myDir}.txt
#rm omega_list_LZinduced${induced_cutoff}.ensemblOrthos.i${myDir}.txt
#cat /mnt/beegfs/ek112884/mus_expression_analysis/${myDir}/INDUCED_CUTOFF${induced_cutoff}/gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm1.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_induced${induced_cutoff}_no_paml.ensemblOrthos.${myDir}.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZinduced${induced_cutoff}.ensemblOrthos.${myDir}.txt; done
#
#rm RS_induced${induced_cutoff}_no_paml.ensemblOrthos.${myDir}.txt
#rm omega_list_RSinduced${induced_cutoff}.ensemblOrthos.${myDir}.txt
#cat /mnt/beegfs/ek112884/mus_expression_analysis/${myDir}/INDUCED_CUTOFF${induced_cutoff}/gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm1.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_induced${induced_cutoff}_no_paml.ensemblOrthos.${myDir}.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RSinduced${induced_cutoff}.ensemblOrthos.${myDir}.txt; done

echo "Done!"

###OLD###

#All commands append to text files, so make sure they aren't already in directory
#rm omega_list_RS.ensemblOrthos.txt
#rm RS_expressed_no_paml.ensemblOrthos.txt
#rm omega_list_RSinduced.ensemblOrthos.txt
#rm RS_induced_no_paml.ensemblOrthos.txt
#rm omega_list_LZ.ensemblOrthos.txt
#rm LZ_expressed_no_paml.ensemblOrthos.txt
#rm omega_list_LZinduced.ensemblOrthos.txt
#rm LZ_induced_no_paml.ensemblOrthos.txt
#rm omega_list.ensemblOrthos.txt
#rm no_paml.ensemblOrthos.txt
#
#cat ../../mus_expression_analysis/MULTI_MAP/gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_expressed_no_paml.ensemblOrthos.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RS.ensemblOrthos.txt; done
#
#cat ../../mus_expression_analysis/MULTI_MAP/gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_induced_no_paml.ensemblOrthos.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RSinduced.ensemblOrthos.txt; done
#
#cat ../../mus_expression_analysis/MULTI_MAP/gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_expressed_no_paml.ensemblOrthos.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZ.ensemblOrthos.txt; done
#
#cat ../../mus_expression_analysis/MULTI_MAP/gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_induced_no_paml.ensemblOrthos.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZinduced.ensemblOrthos.txt; done
#
#tail -n +2 ../../mus_expression_analysis/orthologs/ensembl_orthos_pairwiseOne2one_combined.txt | awk '{print $1}' | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> no_paml.ensemblOrthos.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list.ensemblOrthos.txt; done
#
##ls -d PAML_RESULTS/CODEML_OUTPUT_ENSMUSG00000* | while read dir; do gene=$(echo "${dir}" | cut -d "_" -f 4); ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list.ensemblOrthos.txt; done
