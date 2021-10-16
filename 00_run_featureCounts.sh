#!/bin/bash
#PURPOSE: run featureCounts for all tophat output bams (RNAseq cell sort mapping output) to get raw counts for exons
#
# Job name:
#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=4 # Number of cores per MPI rank (ie number of threads, I think)
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
map_path="/mnt/beegfs/ek112884/mus_expression_analysis/"
annot_path="/mnt/beegfs/ek112884/mus_expression_analysis/ref/"
out_path="/mnt/beegfs/ek112884/mus_expression_analysis/WHOLE_GENOME_NO_MULTIMAP/"

#M.m.musculus
#PWK
#ls ${map_path}PPPP*tophat_out/accepted_hits.bam | while read file; do
#	name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#	echo "${name}"
#	#WITH multi-mapping
#	#featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##CZII
#ls ${map_path}CCCC*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#        #WITH multi-mapping
#	#featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##MBS
#ls ${map_path}MBS*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_pwkphj.PWK_PhJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##M.m.domesticus
##BIK
#ls ${map_path}BIK*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##DGA
#ls ${map_path}DGA*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##LEWES
#ls ${map_path}LLLL*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##WSB
#ls ${map_path}WWWW*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_musculus_wsbeij.WSB_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done
#
##M.spretus
#ls ${map_path}S*tophat_out/accepted_hits.bam | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
#        echo "${name}"
#	#WITH multi-mapping
#        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_spretus.SPRET_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#	#WITHOUT multi-mapping
#	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_spretus.SPRET_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
#done

#M.pahari
ls ${map_path}PAH*tophat_out/accepted_hits.bam | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 6 | cut -d "_" -f 1)
        echo "${name}"
	#WITH multi-mapping
        #featureCounts -T 4 -M -t exon -g gene_id -a ${annot_path}Mus_pahari.PAHARI_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
	#WITHOUT multi-mapping
	featureCounts -T 4 -t exon -g gene_id -a ${annot_path}Mus_pahari.PAHARI_EiJ_v1.94.gtf.gz -o "${out_path}${name}_counts.txt" "${file}"
done

echo "Done!"
