#!/bin/bash

#PURPOSE: generate alignment and tree files for each gene for codeml input
#If multiple isoforms, align all and then filter for mult of 3; if multiple alignments are mult 3 then take longest
#rm genes_not_mult3.txt
ID=""
#zcat ../cds_in_longest_txpts.bed.gz | head -28 | while read gene
#zcat ../cds_in_longest_txpts.bed.gz | while read gene
head -11000 genes_for_alignment.txt | tail -1000 | while read gene
do
        temp_ID=${gene}
        #if [[ -d "DONE/CODEML_OUTPUT_${temp_ID}" ]]
        #then
        #        continue
        #fi

        if [[ ${ID} = ${temp_ID} ]]
        then
                continue
        else
                 ID=${temp_ID}
        fi

	echo "${ID}"

	grep "${ID}" PWK_PhJ.cds_for_paml_alignment.bed > temp.bed
        #Get continuous, nonoverlapping CDS for gene
        total_lines=$(wc -l temp.bed | awk '{print $1}')
        l_n=0
        dir=$(awk '{print $6}' temp.bed | head -1)
        longest_iso=0
	echo "dir: ${dir}, total lines: ${total_lines}, LN: ${l_n}"
        #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
        while [[ ${l_n} -lt ${total_lines} ]]; do
        	l_n=$((${l_n}+1))
                head -${l_n} temp.bed | tail -1 > new_temp.bed
                if [ "${dir}" == "-" ]; then
                	lastline=$(awk '{print $2}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                        	#echo "${l}"
                                l_array=($l)
                                if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
                                	#echo "${l_array[2]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[1]}
                                        l_n=$(($l_n + 1))
                               fi
                         done < temp2.bed
                fi
		if [ "${dir}" == "+" ]; then
                	lastline=$(awk '{print $3}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                        	l_array=($l)
                                if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
                                	#echo "${l_array[1]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[2]}
                                        l_n=$(($l_n + 1))
                                fi
                        done < temp2.bed
                fi
	iso_len=$(awk '{print $3-$2}' new_temp.bed | paste -sd+ | bc)
        echo "${iso_len}"
        if [[ ${iso_len} -gt ${longest_iso} ]]; then
                longest_iso=${iso_len}
                cat new_temp.bed > PWK_temp.bed
        fi
	done
	grep "${ID}" WSB_EiJ.cds_for_paml_alignment.bed > temp.bed
        #Get continuous, nonoverlapping CDS for gene
        total_lines=$(wc -l temp.bed | awk '{print $1}')
        l_n=0
        dir=$(awk '{print $6}' temp.bed | head -1)
        longest_iso=0
	echo "dir: ${dir}, total lines: ${total_lines}, LN: ${l_n}"
        #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
        while [[ ${l_n} -lt ${total_lines} ]]; do
                l_n=$((${l_n}+1))
                head -${l_n} temp.bed | tail -1 > new_temp.bed
                if [ "${dir}" == "-" ]; then
                        lastline=$(awk '{print $2}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                #echo "${l}"
                                l_array=($l)
                                if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
                                        #echo "${l_array[2]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[1]}
                                        l_n=$(($l_n + 1))
                               fi
                         done < temp2.bed
                fi
		if [ "${dir}" == "+" ]; then
                        lastline=$(awk '{print $3}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                l_array=($l)
                                if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
                                        #echo "${l_array[1]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[2]}
                                        l_n=$(($l_n + 1))
                                fi
                        done < temp2.bed
                fi
	iso_len=$(awk '{print $3-$2}' new_temp.bed | paste -sd+ | bc)
        echo "${iso_len}"
        if [[ ${iso_len} -gt ${longest_iso} ]]; then
                longest_iso=${iso_len}
                cat new_temp.bed > WSB_temp.bed
        fi
        done
	grep "${ID}" SPRET_EiJ.cds_for_paml_alignment.bed > temp.bed
        #Get continuous, nonoverlapping CDS for gene
        total_lines=$(wc -l temp.bed | awk '{print $1}')
        l_n=0
        dir=$(awk '{print $6}' temp.bed | head -1)
	longest_iso=0
        echo "dir: ${dir}, total lines: ${total_lines}, LN: ${l_n}"
        #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
        while [[ ${l_n} -lt ${total_lines} ]]; do
                l_n=$((${l_n}+1))
                head -${l_n} temp.bed | tail -1 > new_temp.bed
                if [ "${dir}" == "-" ]; then
                        lastline=$(awk '{print $2}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                #echo "${l}"
                                l_array=($l)
                                if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
                                        #echo "${l_array[2]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[1]}
                                        l_n=$(($l_n + 1))
                               fi
                         done < temp2.bed
                fi
		if [ "${dir}" == "+" ]; then
                        lastline=$(awk '{print $3}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                l_array=($l)
                                if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
                                        #echo "${l_array[1]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[2]}
                                        l_n=$(($l_n + 1))
                                fi
                        done < temp2.bed
                fi
        iso_len=$(awk '{print $3-$2}' new_temp.bed | paste -sd+ | bc)
	echo "${iso_len}"
	if [[ ${iso_len} -gt ${longest_iso} ]]; then
		longest_iso=${iso_len}
		cat new_temp.bed > SPRET_temp.bed
	fi
        done
	grep "${ID}" PAHARI_EiJ.cds_for_paml_alignment.bed > temp.bed
        #Get continuous, nonoverlapping CDS for gene
        total_lines=$(wc -l temp.bed | awk '{print $1}')
        l_n=0
        dir=$(awk '{print $6}' temp.bed | head -1)
	longest_iso=0
        echo "dir: ${dir}, total lines: ${total_lines}, LN: ${l_n}"
        #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
        while [[ ${l_n} -lt ${total_lines} ]]; do
                l_n=$((${l_n}+1))
                head -${l_n} temp.bed | tail -1 > new_temp.bed
                if [ "${dir}" == "-" ]; then
                        lastline=$(awk '{print $2}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                #echo "${l}"
                                l_array=($l)
                                if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
                                        #echo "${l_array[2]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[1]}
                                        l_n=$(($l_n + 1))
                               fi
                         done < temp2.bed
                fi
		if [ "${dir}" == "+" ]; then
                        lastline=$(awk '{print $3}' temp.bed | head -1)
                        tail -n +${l_n} temp.bed > temp2.bed
                        while read -r l; do
                                l_array=($l)
                                if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
                                        #echo "${l_array[1]} ${lastline}"
                                        echo "${l}" >> new_temp.bed
                                        lastline=${l_array[2]}
                                        l_n=$(($l_n + 1))
                                fi
                        done < temp2.bed
                fi
        iso_len=$(awk '{print $3-$2}' new_temp.bed | paste -sd+ | bc)
        echo "${iso_len}"
        if [[ ${iso_len} -gt ${longest_iso} ]]; then
                longest_iso=${iso_len}
                cat new_temp.bed > PAHARI_temp.bed
        fi
	done
        
#	echo "${ID}"
#        grep "${ID}" PAHARI_EiJ.cds_in_longest_txpts.bed > pahari_temp.bed
#	grep "${ID}" PWK_PhJ.cds_in_longest_txpts.bed > PWK_temp.bed
#	grep "${ID}" SPRET_EiJ.cds_in_longest_txpts.bed > SPRET_temp.bed
#	grep "${ID}" WSB_EiJ.cds_in_longest_txpts.bed > WSB_temp.bed

	file=$(echo "Pahari_EiJ.chromosomes.unplaced.gt2k.fa")
	#echo "${file}"
	name=$(echo "${file}" | cut -d "." -f 1)
	bedtools getfasta -fi "${file}" -bed PAHARI_temp.bed -fo "${name}.temp.fa" -s
	sed -i '1!{/^>/d;}' "${name}.temp.fa"
	sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
	sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
	cat "${name}.temp.fa" >> temp.fa
	rm "${name}.temp.fa"

	file=$(echo "PWK_PhJ.chromosomes.unplaced.gt2k.fa")
        name=$(echo "${file}" | cut -d "." -f 1)
        bedtools getfasta -fi "${file}" -bed PWK_temp.bed -fo "${name}.temp.fa" -s
        sed -i '1!{/^>/d;}' "${name}.temp.fa"
        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
        cat "${name}.temp.fa" >> temp.fa
        rm "${name}.temp.fa"

	file=$(echo "SPRET_EiJ.chromosomes.unplaced.gt2k.fa")
        name=$(echo "${file}" | cut -d "." -f 1)
        bedtools getfasta -fi "${file}" -bed SPRET_temp.bed -fo "${name}.temp.fa" -s
        sed -i '1!{/^>/d;}' "${name}.temp.fa"
        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
        cat "${name}.temp.fa" >> temp.fa
        rm "${name}.temp.fa"

	file=$(echo "WSB_EiJ.chromosomes.unplaced.gt2k.fa")
        name=$(echo "${file}" | cut -d "." -f 1)
        bedtools getfasta -fi "${file}" -bed WSB_temp.bed -fo "${name}.temp.fa" -s
        sed -i '1!{/^>/d;}' "${name}.temp.fa"
        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
        cat "${name}.temp.fa" >> temp.fa
        rm "${name}.temp.fa"

        mkdir CODEML_OUTPUT_${ID}
        cd CODEML_OUTPUT_${ID}
        mv ../temp.fa .

        echo "##########Mafft alignment########################"
        mafft --auto --thread 8 temp.fa > mafft_aligned.fa
        rm temp.fa

        echo "##########AMAS#################################"
        AMAS convert -d dna -f fasta -i mafft_aligned.fa -u phylip -c 8
        rm mafft_aligned.fa

#       echo "#########Append N's so all sequences multiples of 3############"
        #seq_len=$(awk '{print length($0)}' temp.fa | tail -1)
#       seq_len=$(head -1 mafft_aligned.fa-out.phy | awk '{print $2}')
#       echo "${seq_len}"
#       remainder=$((seq_len % 3))
 #       echo "${remainder}"
  #      if [ ${remainder} -eq 1 ]
   #     then
    #            sed -i '/^[A-Za-z]/s/$/NN/' mafft_aligned.fa-out.phy
                #sed -i '/^[A-Za-z]/s/^/NN/' temp.fa
#                sed -i "/^[0-9]/s/[ ][0-9]*/ $((${seq_len}+2))/" mafft_aligned.fa-out.phy
#        fi
#        if [ ${remainder} -eq 2 ]
#        then
#                sed -i '/^[A-Za-z]/s/$/N/' mafft_aligned.fa-out.phy
                #sed -i '/^[A-Za-z]/s/^/N/' temp.fa
#                sed -i "/^[0-9]/s/[ ][0-9]*/ $((${seq_len}+1))/" mafft_aligned.fa-out.phy
#        fi


        mv mafft_aligned.fa-out.phy final_alignment.phy
        #rm mafft_aligned.fa-out.phy
	
	echo "alignment length:"
	head -1 final_alignment.phy | awk '{print $2}'
	rem=$(head -1 final_alignment.phy | awk '{print $2%3}')
	if [[ ${rem} -eq 0 ]]; then
	        echo "##########IQ tree###############################"
	        sta=$(echo "${gene}" | awk '{print $2}')
	        sto=$(echo "${gene}" | awk '{print $3}')
		if [[ $sta > 1 && $sto > 1 ]]
	        then
	                length=$((${sto}-${sta}))
	                div=$((${length}/100000+1))
	                iqtree-omp -s final_alignment.phy -nt ${div} -redo -mredo
	        else
	                iqtree-omp -s final_alignment.phy -nt 2 -redo -mredo
	        fi
		cd /scratch/1/kopania/PAML_runs/TEMP1_WG
	else
		cd /scratch/1/kopania/PAML_runs/TEMP1_WG
		rm -r CODEML_OUTPUT_${ID}
		echo "${ID}" >> genes_not_mult3.txt
	fi

        cd /scratch/1/kopania/PAML_runs/TEMP1_WG

done
 
