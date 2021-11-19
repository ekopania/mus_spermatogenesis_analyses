#!/bin/bash

#PURPOSE: Lots of genes have alignments not multiples of 3 after aligning all four species from longest txpts; go through these genes and see if other transcripts (not necessarily longest) have alignments that are mult 3
#Note: doing this in a hacky where where I manually check the first CDS, then the second, etc

cat genes_list_not_mult3.txt | while read gene
do
        temp_ID=${gene}
        if [[ -d "PAML_RESULTS/CODEML_OUTPUT_${temp_ID}" ]]
        then
        	continue
        fi

        if [[ ${ID} = ${temp_ID} ]]
        then
                continue
        else
                 ID=${temp_ID}
        fi
        echo "${ID}"
        grep "${ID}" PAHARI_EiJ.cds_for_paml_alignment.uniq.bed | head -1 > pahari_temp.bed
        grep "${ID}" PWK_PhJ.cds_for_paml_alignment.uniq.bed | head -1 > PWK_temp.bed
        grep "${ID}" SPRET_EiJ.cds_for_paml_alignment.uniq.bed | head -1 > SPRET_temp.bed
        grep "${ID}" WSB_EiJ.cds_for_paml_alignment.uniq.bed | head -1 > WSB_temp.bed

        file=$(echo "Pahari_EiJ.chromosomes.unplaced.gt2k.fa")
        #echo "${file}"
        name=$(echo "${file}" | cut -d "." -f 1)
        bedtools getfasta -fi "${file}" -bed pahari_temp.bed -fo "${name}.temp.fa" -s
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

	mv mafft_aligned.fa-out.phy final_alignment.phy
        rm mafft_aligned.fa-out.phy

	r=$(head -1 final_alignment.phy | awk '{print ($2 % 3)}')
	if [ $r -eq 0 ]
	then
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

        	cd /scratch/1/kopania/PAML_runs/WHOLE_GENOMES
	else
		cd /scratch/1/kopania/PAML_runs/WHOLE_GENOMES
		rm -r CODEML_OUTPUT_${ID}
	fi
done

