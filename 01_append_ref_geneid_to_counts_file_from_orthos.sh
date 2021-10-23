#!/bin/bash
#
#PURPOSE: all of the *counts.txt files output from feature counts have species-specific gene IDs, but they need to be matched to the same gene ID for EdgeR to process them, so appending the gene ID of the mouse reference ortholog to all of these
#
# Job name:
#SBATCH --job-name=append_orthoID
#SBATCH --output=append_orthoID-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-3 # run on a specific node
#
## Command(s) to run:
ortho_file="/mnt/beegfs/ek112884/mus_expression_analysis/orthologs/ensembl_orthos_pairwiseOne2one_combined.txt"

#Multi map, count reads for each mapping (featurecounts -M)
#file_path="/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/"
#Multi map, fractional
file_path="/mnt/beegfs/ek112884/mus_expression_analysis/WHOLE_GENOME_MULTIMAP_FRACTIONAL/"
#Do not count multimapping reads
#file_path="/mnt/beegfs/ek112884/mus_expression_analysis/WHOLE_GENOME_NO_MULTIMAP/"

echo "Working on files in this directory: ${file_path}"

ls ${file_path}BIK*counts.txt | while read file; do
	tail -n +3 "${file}" > temp1.txt
	join -1 2 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 2 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
	name=$(echo "${file}" | cut -d "." -f 1)
	echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
	cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}DGA*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
        join -1 2 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 2 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
	name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}LLLL*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 2 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 2 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
        echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}WWWW*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
        join -1 2 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 2 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
	name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}CCCC*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 3 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 3 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}MBS*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 3 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 3 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt	
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"        
done

ls ${file_path}PPPP*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 3 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 3 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt	
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"        
done

ls ${file_path}S*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 4 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 4 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

ls ${file_path}PAHNew*counts.txt | while read file; do
        tail -n +3 "${file}" > temp1.txt
	join -1 5 -2 1 -o 0,1.1,2.2,2.3,2.4,2.5,2.6,2.7 -t $'\t' <(sort -k 5 ${ortho_file}) <(sort -k 1 temp1.txt) > temp2.txt
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
	echo -e "Geneid_spec\tGeneid\tChr\tStart\tEnd\tStrand\tLength\tCounts" > "${name}.orthoIDappended.txt"
        cat temp2.txt >> "${name}.orthoIDappended.txt"
done

rm temp1.txt
rm temp2.txt

echo "Done!"
