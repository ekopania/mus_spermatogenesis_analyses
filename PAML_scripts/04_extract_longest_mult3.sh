#!/bin/bash

#PURPOSE: For genes with multiple alignments that are multiples of 3, take the longest and delete directories containing other alignments

mkdir NOT_LONGEST
mkdir MORE_GAPS
mkdir REPLICATES
#Loop through all genes (all CODEML directories)
ls -d *_1 | while read file; do 
	dirname=$(echo "${file}" | cut -d "_" -f 1-3)
	echo "${dirname}"
	longest=0
	least_gaps=9223372036854775807
	#Loop through all alignments for a given gene
	while read -r dir; do
		#If longest alignment, keep here and inc $longest, otherwise move to "NOT_LONGEST" dir
		head -1 "${dir}/final_alignment.phy"
		align_len=$(head -1 "${dir}/final_alignment.phy" | awk '{print $2}')
		if [ $align_len -ge $longest ]; then
			longest=${align_len}
			echo "Longest: ${longest}"
		else
			mv "${dir}" NOT_LONGEST
		fi
	done < <(ls -d ${dirname}_*)
	#This approach will always save the first alignment of a gene as a "longest alignment"; check to make sure it actually is a longest alignment
        first_align_len=$(head -1 "${dirname}_1/final_alignment.phy" | awk '{print $2}')
        if [ $first_align_len -lt $longest ]; then
        	mv ${dirname}_1 NOT_LONGEST
        fi
	#Loop through remaining (longest) alignments for a given gene and save the one w/ fewest gaps
	while read -r dir; do
		ngaps=$(echo $(cat "${dir}/final_alignment.phy" | wc -c) - $(cat "${dir}/final_alignment.phy" | tr -d "-" | wc -c) | bc)
                echo "${ngaps}"
                #If least amount of gaps so far, keep and dec $least_gaps, otherwise move to "MORE_GAPS" dir
                if [ $ngaps -le $least_gaps ]; then
                	least_gaps=${ngaps}
                        echo "Least gaps: ${least_gaps}"
		else
			mv "${dir}" MORE_GAPS
                fi
	done < <(ls -d ${dirname}_*)
	#This approach will always save the first alignment of a gene as having the fewest gaps; check to make sure it actually has fewest gaps
	d=$(find "./" -maxdepth 1 -type d -name "${dirname}_*" | head -1)
	echo "d: ${d}"
	first_align_ngaps=$(echo $(cat "${d}/final_aligment.phy" | wc -c) - $(cat "${d}/final_alignment.phy" | tr -d "-" | wc -c) | bc)
	echo "First align ngaps: ${first_align_ngaps}"
	if [ $first_align_ngaps -lt $least_gaps ]; then
		mv "${d}" MORE_GAPS
	fi
	#May still be multiple alignments left; make sure they're exactly the same and then keep first one
	num_aligns=$(ls -d ${dirname}_* | wc -l)
	check_dif=$(diff ${dirname}_*/final_alignment.phy)
	if [ $num_aligns -gt 1 ] && [ ! $check_dif ]; then
		ls -d ${dirname}* | tail -$(($num_aligns-1)) | while read rep_dir; do
			mv "${rep_dir}" REPLICATES
		done
	fi
done
