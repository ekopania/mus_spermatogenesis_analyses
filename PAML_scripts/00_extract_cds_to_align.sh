#!/bin/bash

#PWK
#rm PWK_PhJ.cds_for_paml_alignment.bed
##head -4 PWK_uniq_geneids.txt | while read ID; do
#cat PWK_uniq_geneids.txt | while read ID; do
#	rm new_temp.bed
#	rm new.fa
#	#echo "${line}"
#	#line_array=($line)
#	#ID=$(echo "${line_array[4]}")
#	echo "${ID}"
#	#Check if in 1:1 orthos
#	if grep -q "${ID}" /scratch/1/kopania/mus_expression_analysis/orthologs/final_combined_orthos.txt; then
#		#echo "${ID}"
#		grep "${ID}" /scratch/1/kopania/mus_expression_analysis/ref/PWK_PhJ.bed > temp.bed
#		#Get continuous, nonoverlapping CDS for gene
#		total_lines=$(wc -l temp.bed | awk '{print $1}')
#		l_n=0
#		dir=$(awk '{print $6}' temp.bed | head -1)
#		#echo "dir: ${dir}, total lines: ${total_lines}, LN: ${l_n}"
#		#Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
#		while [[ ${l_n} -lt ${total_lines} ]]; do
#			l_n=$((${l_n}+1))
#			head -${l_n} temp.bed | tail -1 > new_temp.bed
#			if [ "${dir}" == "-" ]; then
#				lastline=$(awk '{print $2}' temp.bed | head -1)
#				tail -n +${l_n} temp.bed > temp2.bed
#				while read -r l; do
#					#echo "${l}"
#					l_array=($l)
#					if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
#						#echo "${l_array[2]} ${lastline}"
#						echo "${l}" >> new_temp.bed
#						lastline=${l_array[1]}
#						l_n=$(($l_n + 1))
#					fi
#				done < temp2.bed
#			fi
#			if [ "${dir}" == "+" ]; then
#				lastline=$(awk '{print $3}' temp.bed | head -1)
#				tail -n +${l_n} temp.bed > temp2.bed
#				while read -r l; do
#					l_array=($l)
#					if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
#						#echo "${l_array[1]} ${lastline}"
#						echo "${l}" >> new_temp.bed
#						lastline=${l_array[2]}
#						l_n=$(($l_n + 1))
#					fi
#				done < temp2.bed
#			fi
#			#Use CDS coordinates from bedfile to get actual sequences from fasta and translate to aa seq to check for missing start codons or early stop codons
#			bedtools getfasta -fi /scratch/1/kopania/mus_expression_analysis/ref/PWK_PhJ.chromosomes.unplaced.gt2k.fa -bed new_temp.bed -fo temp.fa -s
#			echo ">PWK" > new.fa
#			grep -v ">" temp.fa | tr '\n' ' ' | sed -e 's/ //g' >> new.fa		
#	                AMAS translate -f phylip -d dna -i new.fa --out-format phylip > translation.out
#	                aa_seq=$(tail -1 translated_new.fa-out.phy-out.phy | awk '{print $2}')
#			#echo "${aa_seq}"
#	                #Check for no early stop codons and start codon exists
#	                if ! grep -q "stop" translation.out && [[ $aa_seq =~ ^M ]]; then
#	                        cat new_temp.bed >> PWK_PhJ.cds_for_paml_alignment.bed
#	                fi
#		done
#	fi
#done

#WSB
#rm WSB_EiJ.cds_for_paml_alignment.bed
#cat WSB_uniq_geneids.txt | while read ID; do
#        rm new_temp.bed
#        rm new.fa
#        echo "${ID}"
#        #Check if in 1:1 orthos
#        if grep -q "${ID}" /scratch/1/kopania/mus_expression_analysis/orthologs/final_combined_orthos.txt; then
#                grep "${ID}" /scratch/1/kopania/mus_expression_analysis/ref/WSB_EiJ.bed > temp.bed
#                #Get continuous, nonoverlapping CDS for gene
#                total_lines=$(wc -l temp.bed | awk '{print $1}')
#                l_n=0
#                dir=$(awk '{print $6}' temp.bed | head -1)
#                #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
#                while [[ ${l_n} -lt ${total_lines} ]]; do
#                        l_n=$((${l_n}+1))
#                        head -${l_n} temp.bed | tail -1 > new_temp.bed
#                        if [ "${dir}" == "-" ]; then
#                                lastline=$(awk '{print $2}' temp.bed | head -1)
#                                tail -n +${l_n} temp.bed > temp2.bed
#                                while read -r l; do
#					l_array=($l)
#                                        if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
#                                                echo "${l}" >> new_temp.bed
#                                                lastline=${l_array[1]}
#                                                l_n=$(($l_n + 1))
#                                        fi
#                                done < temp2.bed
#                        fi
#                        if [ "${dir}" == "+" ]; then
#                                lastline=$(awk '{print $3}' temp.bed | head -1)
#                                tail -n +${l_n} temp.bed > temp2.bed
#                                while read -r l; do
#                                        l_array=($l)
#                                        if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
#                                                echo "${l}" >> new_temp.bed
#                                                lastline=${l_array[2]}
#                                                l_n=$(($l_n + 1))
#                                        fi
#                                done < temp2.bed
#                        fi
#                        #Use CDS coordinates from bedfile to get actual sequences from fasta and translate to aa seq to check for missing start codons or early stop codons
#                        bedtools getfasta -fi /scratch/1/kopania/mus_expression_analysis/ref/WSB_EiJ.chromosomes.unplaced.gt2k.fa -bed new_temp.bed -fo temp.fa -s
#                        echo ">WSB" > new.fa
#                        grep -v ">" temp.fa | tr '\n' ' ' | sed -e 's/ //g' >> new.fa
#                        AMAS translate -f phylip -d dna -i new.fa --out-format phylip > translation.out
#                        aa_seq=$(tail -1 translated_new.fa-out.phy-out.phy | awk '{print $2}')
#			#Check for no early stop codons and start codon exists
#                        if ! grep -q "stop" translation.out && [[ $aa_seq =~ ^M ]]; then
#                                cat new_temp.bed >> WSB_EiJ.cds_for_paml_alignment.bed
#                        fi
#                done
#        fi
#done

#SPRET
#rm SPRET_EiJ.cds_for_paml_alignment.bed
#cat SPRET_uniq_geneids.txt | while read ID; do
#        rm new_temp.bed
#        rm new.fa
#        echo "${ID}"
#        #Check if in 1:1 orthos
#        if grep -q "${ID}" /scratch/1/kopania/mus_expression_analysis/orthologs/final_combined_orthos.txt; then
#                grep "${ID}" /scratch/1/kopania/mus_expression_analysis/ref/SPRET_EiJ.bed > temp.bed
#                #Get continuous, nonoverlapping CDS for gene
#                total_lines=$(wc -l temp.bed | awk '{print $1}')
#                l_n=0
#                dir=$(awk '{print $6}' temp.bed | head -1)
#                #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
#                while [[ ${l_n} -lt ${total_lines} ]]; do
#                        l_n=$((${l_n}+1))
#                        head -${l_n} temp.bed | tail -1 > new_temp.bed
#                        if [ "${dir}" == "-" ]; then
#                                lastline=$(awk '{print $2}' temp.bed | head -1)
#                                tail -n +${l_n} temp.bed > temp2.bed
#                                while read -r l; do
#                                        l_array=($l)
#                                        if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
#                                                echo "${l}" >> new_temp.bed
#                                                lastline=${l_array[1]}
#                                                l_n=$(($l_n + 1))
#                                        fi
#                                done < temp2.bed
#                        fi
#			if [ "${dir}" == "+" ]; then
#                                lastline=$(awk '{print $3}' temp.bed | head -1)
#                                tail -n +${l_n} temp.bed > temp2.bed
#                                while read -r l; do
#                                        l_array=($l)
#                                        if [[ "${l_array[1]}" -gt "${lastline}" ]]; then
#                                                echo "${l}" >> new_temp.bed
#                                                lastline=${l_array[2]}
#                                                l_n=$(($l_n + 1))
#                                        fi
#                                done < temp2.bed
#                        fi
#                        #Use CDS coordinates from bedfile to get actual sequences from fasta and translate to aa seq to check for missing start codons or early stop codons
#                        bedtools getfasta -fi /scratch/1/kopania/mus_expression_analysis/ref/SPRET_EiJ.chromosomes.unplaced.gt2k.fa -bed new_temp.bed -fo temp.fa -s
#                        echo ">SPRET" > new.fa
#                        grep -v ">" temp.fa | tr '\n' ' ' | sed -e 's/ //g' >> new.fa
#                        AMAS translate -f phylip -d dna -i new.fa --out-format phylip > translation.out
#                        aa_seq=$(tail -1 translated_new.fa-out.phy-out.phy | awk '{print $2}')
#                        #Check for no early stop codons and start codon exists
#                        if ! grep -q "stop" translation.out && [[ $aa_seq =~ ^M ]]; then
#                                cat new_temp.bed >> SPRET_EiJ.cds_for_paml_alignment.bed
#                        fi
#                done
#        fi
#done

#PAHARI
rm PAHARI_EiJ.cds_for_paml_alignment.bed
cat PAHARI_uniq_geneids.txt | while read ID; do
        rm new_temp.bed
        rm new.fa
        echo "${ID}"
        #Check if in 1:1 orthos
        if grep -q "${ID}" /scratch/1/kopania/mus_expression_analysis/orthologs/final_combined_orthos.txt; then
                grep "${ID}" /scratch/1/kopania/mus_expression_analysis/ref/PAHARI_EiJ.bed > temp.bed
                #Get continuous, nonoverlapping CDS for gene
                total_lines=$(wc -l temp.bed | awk '{print $1}')
                l_n=0
                dir=$(awk '{print $6}' temp.bed | head -1)
                #Loop through each transcript/CDS (start over when coordinates get higher or lower dependiing on direction (- or + respectively), do this for each isoform of gene
                while [[ ${l_n} -lt ${total_lines} ]]; do
                        l_n=$((${l_n}+1))
                        head -${l_n} temp.bed | tail -1 > new_temp.bed
                        if [ "${dir}" == "-" ]; then
                                lastline=$(awk '{print $2}' temp.bed | head -1)
                                tail -n +${l_n} temp.bed > temp2.bed
                                while read -r l; do
                                        l_array=($l)
                                        if [[ "${l_array[2]}" -lt "${lastline}" ]]; then
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
                                                echo "${l}" >> new_temp.bed
                                                lastline=${l_array[2]}
                                                l_n=$(($l_n + 1))
                                        fi
                                done < temp2.bed
                        fi
                        #Use CDS coordinates from bedfile to get actual sequences from fasta and translate to aa seq to check for missing start codons or early stop codons
                        bedtools getfasta -fi /scratch/1/kopania/mus_expression_analysis/ref/Pahari_EiJ.chromosomes.unplaced.gt2k.fa -bed new_temp.bed -fo temp.fa -s
                        echo ">PAHARI" > new.fa
                        grep -v ">" temp.fa | tr '\n' ' ' | sed -e 's/ //g' >> new.fa
                        AMAS translate -f phylip -d dna -i new.fa --out-format phylip > translation.out
                        aa_seq=$(tail -1 translated_new.fa-out.phy-out.phy | awk '{print $2}')
                        #Check for no early stop codons and start codon exists
                        if ! grep -q "stop" translation.out && [[ $aa_seq =~ ^M ]]; then
                                cat new_temp.bed >> PAHARI_EiJ.cds_for_paml_alignment.bed
                        fi
                done
        fi
done
