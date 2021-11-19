#!/bin/bash

#PURPOSE: Loop through genes not included in original alignments (genes_not_mult3.txt) and see if there are alternative isoforms that can be included

#Check if gene included in all four species bedfiles (if not, go on to next gene)
rm genes_not_in_all_four.txt
cat genes_not_mult3.txt | while read ID; do
	echo "${ID}"
	if (grep -q "${ID}" PWK_PhJ.cds_for_paml_alignment.bed) && (grep -q "${ID}" WSB_EiJ.cds_for_paml_alignment.bed) && (grep -q "${ID}" SPRET_EiJ.cds_for_paml_alignment.bed) && (grep -q "${ID}" PAHARI_EiJ.cds_for_paml_alignment.bed); then 
		grep "${ID}" PWK_PhJ.cds_for_paml_alignment.bed > temp.bed
	        #Get continuous, nonoverlapping CDS for gene
	        total_lines=$(wc -l temp.bed | awk '{print $1}')
	        l_n=0
	        dir=$(awk '{print $6}' temp.bed | head -1)
	        longest_iso=0
		PWK_num_isos=0
		rm PWK_temp*.bed
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
		PWK_num_isos=$(($PWK_num_isos + 1))
		cat new_temp.bed > PWK_temp_${PWK_num_isos}.bed
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
		WSB_num_isos=0
		rm WSB_temp*.bed
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
		WSB_num_isos=$(($WSB_num_isos + 1))
		cat new_temp.bed > WSB_temp_${WSB_num_isos}.bed
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
		SPRET_num_isos=0
		rm SPRET_temp*.bed
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
		SPRET_num_isos=$(($SPRET_num_isos + 1))
		cat new_temp.bed > SPRET_temp_${SPRET_num_isos}.bed
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
		PAHARI_num_isos=0
		rm PAHARI_temp*.bed
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
	        PAHARI_num_isos=$(($PAHARI_num_isos + 1))
		cat new_temp.bed > PAHARI_temp_${PAHARI_num_isos}.bed
		iso_len=$(awk '{print $3-$2}' new_temp.bed | paste -sd+ | bc)
	        echo "${iso_len}"
	        if [[ ${iso_len} -gt ${longest_iso} ]]; then
	                longest_iso=${iso_len}
	                cat new_temp.bed > PAHARI_temp.bed
	        fi
	        done
		#Do all possible alignments of all isoforms of genes across the four species, keep any that are mult3
		inc=1
		ls PWK_temp_*.bed | while read a; do
			ls WSB_temp_*.bed | while read b; do
				ls SPRET_temp_*.bed | while read c; do
					ls PAHARI_temp_*.bed | while read d; do
						file=$(echo "Pahari_EiJ.chromosomes.unplaced.gt2k.fa")
					        #echo "${file}"
					        name=$(echo "${file}" | cut -d "." -f 1)
					        bedtools getfasta -fi "${file}" -bed "${d}" -fo "${name}.temp.fa" -s
					        sed -i '1!{/^>/d;}' "${name}.temp.fa"
					        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
					        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
					        cat "${name}.temp.fa" >> temp.fa
					        rm "${name}.temp.fa"
					
					        file=$(echo "PWK_PhJ.chromosomes.unplaced.gt2k.fa")
					        name=$(echo "${file}" | cut -d "." -f 1)
					        bedtools getfasta -fi "${file}" -bed "${a}" -fo "${name}.temp.fa" -s
					        sed -i '1!{/^>/d;}' "${name}.temp.fa"
					        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
					        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
					        cat "${name}.temp.fa" >> temp.fa
					        rm "${name}.temp.fa"
					
					        file=$(echo "SPRET_EiJ.chromosomes.unplaced.gt2k.fa")
					        name=$(echo "${file}" | cut -d "." -f 1)
					        bedtools getfasta -fi "${file}" -bed "${c}" -fo "${name}.temp.fa" -s
					        sed -i '1!{/^>/d;}' "${name}.temp.fa"
					        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
					        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
					        cat "${name}.temp.fa" >> temp.fa
					        rm "${name}.temp.fa"

						file=$(echo "WSB_EiJ.chromosomes.unplaced.gt2k.fa")
					        name=$(echo "${file}" | cut -d "." -f 1)
					        bedtools getfasta -fi "${file}" -bed "${b}" -fo "${name}.temp.fa" -s
					        sed -i '1!{/^>/d;}' "${name}.temp.fa"
					        sed -i ':a;/^[A-Z]/{N;s/\n//;ba}' "${name}.temp.fa"
					        sed -i "s/>.*/>${name}/g" "${name}.temp.fa"
					        cat "${name}.temp.fa" >> temp.fa
					        rm "${name}.temp.fa"
					
					        mkdir CODEML_OUTPUT_${ID}_${inc}
					        cd CODEML_OUTPUT_${ID}_${inc}
					        mv ../temp.fa .
					
					        echo "##########Mafft alignment########################"
					        mafft --auto --thread 8 temp.fa > mafft_aligned.fa
					        rm temp.fa
					
					        echo "##########AMAS#################################"
					        AMAS convert -d dna -f fasta -i mafft_aligned.fa -u phylip -c 8
					        rm mafft_aligned.fa
						
						mv mafft_aligned.fa-out.phy final_alignment.phy
					
					        echo "alignment length:"
					        head -1 final_alignment.phy | awk '{print $2}'
					        rem=$(head -1 final_alignment.phy | awk '{print $2%3}')
					        if [[ ${rem} -eq 0 ]]; then
							inc=$(($inc + 1))
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
					                cd /scratch/1/kopania/PAML_runs/TEMP2_WG
					        else
					                cd /scratch/1/kopania/PAML_runs/TEMP2_WG
					                rm -r CODEML_OUTPUT_${ID}_${inc}
					        fi
					
					        cd /scratch/1/kopania/PAML_runs/TEMP2_WG
					done
				done
			done
		done
		echo "${PWK_num_isos} ${WSB_num_isos} ${SPRET_num_isos} ${PAHARI_num_isos}"
	else
		echo "${ID}" >> genes_not_in_all_four.txt
	fi
done
