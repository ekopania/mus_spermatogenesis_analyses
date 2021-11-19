#!/bin/bash

#PURPOSE: Some of the gene IDs from the mouse reference bed file are not present in the bedfiles for all four species (PWK,WSB,SPRET,Pahari). Could be not annotated or not ortho exists. We only want genes present in all four for dN/dS analyses so this script identifies and separates out genes that do not have an ortho in all four.

mkdir NO_ORTHO_IN_ALL_FOUR 
ls -d CODEML_OUTPUT_* | while read dir; do
	if [ ! -f "${dir}/final_alignment.phy" ]; then
		mv "${dir}" NO_ORTHO_IN_ALL_FOUR
	else  
		nlines=$(awk '{print $1}' <(wc -l "${dir}/final_alignment.phy"))
		echo "${nlines}"
		if [ "${nlines}" -ne 5 ]; then
			mv "${dir}" NO_ORTHO_IN_ALL_FOUR
		fi
	fi
done
