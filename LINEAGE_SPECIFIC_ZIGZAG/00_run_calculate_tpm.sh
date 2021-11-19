#!/bin/bash
#PURPOSE: Wrapper to run 01_calculate_tpm.r using slurm
#
# Job name:
#SBATCH --job-name=calculate_tpm
#SBATCH --output=calculate_tpm_output.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=192000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#Run R script to calculate tpm
#BIK LZ
#ls ../BIK*LZ_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 BIK4665-1M-LZ.tpm.txt) <(sort -k 1 BIK4665-2M-LZ.tpm.txt) > BIK_LZ.combined.tpm.txt
#echo -e "\tBIK4665-1M-LZ.tpm.txt\tBIK4665-2M-LZ.tpm.txt" | cat - BIK_LZ.combined.tpm.txt > BIK_LZ.combined.tpm.header.txt
#rm BIK_LZ.combined.tpm.txt
#mv BIK_LZ.combined.tpm.header.txt BIK_LZ.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff BIK4665-1M-LZ.length.txt BIK4665-2M-LZ.length.txt
#sort BIK4665-1M-LZ.length.txt > BIK_LZ.length.txt
#
##BIK RS
#ls ../BIK*RS_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 BIK4665-1M-RS.tpm.txt) <(sort -k 1 BIK4665-2M-RS.tpm.txt) > BIK_RS.combined.tpm.txt
#echo -e "\tBIK4665-1M-RS.tpm.txt\tBIK4665-2M-RS.tpm.txt" | cat - BIK_RS.combined.tpm.txt > BIK_RS.combined.tpm.header.txt
#rm BIK_RS.combined.tpm.txt
#mv BIK_RS.combined.tpm.header.txt BIK_RS.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff BIK4665-1M-RS.length.txt BIK4665-2M-RS.length.txt
#sort BIK4665-1M-RS.length.txt > BIK_RS.length.txt
#
##DGA LZ
#ls ../DGA*LZ_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 DGA5406-1M-LZ.tpm.txt) <(sort -k 1 DGA5406-2M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 DGA5406-3M-LZ.tpm.txt) > DGA_LZ.combined.tpm.txt
#echo -e "\tDGA5406-1M-LZ.tpm.txt\tDGA5406-2M-LZ.tpm.txt\tDGA5406-3M-LZ.tpm.txt" | cat - DGA_LZ.combined.tpm.txt > DGA_LZ.combined.tpm.header.txt
#rm DGA_LZ.combined.tpm.txt
#mv DGA_LZ.combined.tpm.header.txt DGA_LZ.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff DGA5406-1M-LZ.length.txt DGA5406-2M-LZ.length.txt
#diff DGA5406-1M-LZ.length.txt DGA5406-3M-LZ.length.txt
#sort DGA5406-1M-LZ.length.txt > DGA_LZ.length.txt
#
##DGA RS
#ls ../DGA*RS_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 DGA5406-1M-RS.tpm.txt) <(sort -k 1 DGA5406-2M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 DGA5406-3M-RS.tpm.txt) > DGA_RS.combined.tpm.txt
#echo -e "\tDGA5406-1M-RS.tpm.txt\tDGA5406-2M-RS.tpm.txt\tDGA5406-3M-RS.tpm.txt" | cat - DGA_RS.combined.tpm.txt > DGA_RS.combined.tpm.header.txt
#rm DGA_RS.combined.tpm.txt
#mv DGA_RS.combined.tpm.header.txt DGA_RS.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff DGA5406-1M-RS.length.txt DGA5406-2M-RS.length.txt
#diff DGA5406-1M-RS.length.txt DGA5406-3M-RS.length.txt
#sort DGA5406-1M-RS.length.txt > DGA_RS.length.txt
#
##LLLL LZ
#ls ../LLLL*LZ_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 LLLL125-1M-LZ.tpm.txt) <(sort -k 1 LLLL125-2M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 LLLL125-3M-LZ.tpm.txt) > LLLL_LZ.combined.tpm.txt
#echo -e "\tLLLL125-1M-LZ.tpm.txt\tLLLL125-2M-LZ.tpm.txt\tLLLL125-3M-LZ.tpm.txt" | cat - LLLL_LZ.combined.tpm.txt > LLLL_LZ.combined.tpm.header.txt
#rm LLLL_LZ.combined.tpm.txt
#mv LLLL_LZ.combined.tpm.header.txt LLLL_LZ.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff LLLL125-1M-LZ.length.txt LLLL125-2M-LZ.length.txt
#diff LLLL125-1M-LZ.length.txt LLLL125-3M-LZ.length.txt
#sort LLLL125-1M-LZ.length.txt > LLLL_LZ.length.txt
#
##LLLL RS
#ls ../LLLL*RS_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 LLLL125-1M-RS.tpm.txt) <(sort -k 1 LLLL125-2M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 LLLL125-3M-RS.tpm.txt) > LLLL_RS.combined.tpm.txt
#echo -e "\tLLLL125-1M-RS.tpm.txt\tLLLL125-2M-RS.tpm.txt\tLLLL125-3M-RS.tpm.txt" | cat - LLLL_RS.combined.tpm.txt > LLLL_RS.combined.tpm.header.txt
#rm LLLL_RS.combined.tpm.txt
#mv LLLL_RS.combined.tpm.header.txt LLLL_RS.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff LLLL125-1M-RS.length.txt LLLL125-2M-RS.length.txt
#diff LLLL125-1M-RS.length.txt LLLL125-3M-RS.length.txt
#sort LLLL125-1M-RS.length.txt > LLLL_RS.length.txt
#
##WWWW LZ
#ls ../WWWW8*LZ_counts.orthoIDappended.txt | while read file; do
#	name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#	echo "${name}"
#	Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 WWWW87-2M-LZ.tpm.txt) <(sort -k 1 WWWW87-3M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 WWWW87-4M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 WWWW89-8M-LZ.tpm.txt) > WWWW_LZ.combined.tpm.txt
#echo -e "\tWWWW87-2M-LZ.tpm.txt\tWWWW87-3M-LZ.tpm.txt\tWWWW87-4M-LZ.tpm.txt\tWWWW89-8M-LZ.tpm.txt" | cat - WWWW_LZ.combined.tpm.txt > WWWW_LZ.combined.tpm.header.txt
#rm WWWW_LZ.combined.tpm.txt
#mv WWWW_LZ.combined.tpm.header.txt WWWW_LZ.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff WWWW87-2M-LZ.length.txt WWWW87-3M-LZ.length.txt
#diff WWWW87-2M-LZ.length.txt WWWW87-4M-LZ.length.txt
#diff WWWW87-2M-LZ.length.txt WWWW89-8M-LZ.length.txt
#sort WWWW87-2M-LZ.length.txt > WWWW_LZ.length.txt
#
##WWWW RS
#ls ../WWWW8*RS_counts.orthoIDappended.txt | while read file; do
#        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
#        echo "${name}"
#        Rscript 01_calculate_tpm.r "${file}"
#done
##combine tpm output files for same genotype samples
#join -1 1 -2 1 <(sort -k 1 WWWW87-2M-RS.tpm.txt) <(sort -k 1 WWWW87-3M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 WWWW87-4M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 WWWW89-8M-RS.tpm.txt) > WWWW_RS.combined.tpm.txt
#echo -e "\tWWWW87-2M-RS.tpm.txt\tWWWW87-3M-RS.tpm.txt\tWWWW87-4M-RS.tpm.txt\tWWWW89-8M-RS.tpm.txt" | cat - WWWW_RS.combined.tpm.txt > WWWW_RS.combined.tpm.header.txt
#rm WWWW_RS.combined.tpm.txt
#mv WWWW_RS.combined.tpm.header.txt WWWW_RS.combined.tpm.txt
##check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
#diff WWWW87-2M-RS.length.txt WWWW87-3M-RS.length.txt
#diff WWWW87-2M-RS.length.txt WWWW87-4M-RS.length.txt
#diff WWWW87-2M-RS.length.txt WWWW89-8M-RS.length.txt
#sort WWWW87-2M-RS.length.txt > WWWW_RS.length.txt

#CCCC LZ
ls ../CCCC*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 CCCC153-4M-LZ.tpm.txt) <(sort -k 1 CCCC153-5M-LZ.tpm.txt) > CCCC_LZ.combined.tpm.txt
echo -e "\tCCCC153-4M-LZ.tpm.txt\tCCCC153-5M-LZ.tpm.txt" | cat - CCCC_LZ.combined.tpm.txt > CCCC_LZ.combined.tpm.header.txt
rm CCCC_LZ.combined.tpm.txt
mv CCCC_LZ.combined.tpm.header.txt CCCC_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff CCCC153-4M-LZ.length.txt CCCC153-5M-LZ.length.txt
sort CCCC153-4M-LZ.length.txt > CCCC_LZ.length.txt

#CCCC RS
ls ../CCCC*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 CCCC153-4M-RS.tpm.txt) <(sort -k 1 CCCC153-5M-RS.tpm.txt) > CCCC_RS.combined.tpm.txt
echo -e "\tCCCC153-4M-RS.tpm.txt\tCCCC153-5M-RS.tpm.txt" | cat - CCCC_RS.combined.tpm.txt > CCCC_RS.combined.tpm.header.txt
rm CCCC_RS.combined.tpm.txt
mv CCCC_RS.combined.tpm.header.txt CCCC_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff CCCC153-4M-RS.length.txt CCCC153-5M-RS.length.txt
sort CCCC153-4M-RS.length.txt > CCCC_RS.length.txt

#MBS LZ
ls ../MBS*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 MBS4527-1M-LZ.tpm.txt) <(sort -k 1 MBS4527-2M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 MBS4527-3M-LZ.tpm.txt) > MBS_LZ.combined.tpm.txt
echo -e "\tMBS4527-1M-LZ.tpm.txt\tMBS4527-2M-LZ.tpm.txt\tMBS4527-3M-LZ.tpm.txt" | cat - MBS_LZ.combined.tpm.txt > MBS_LZ.combined.tpm.header.txt
rm MBS_LZ.combined.tpm.txt
mv MBS_LZ.combined.tpm.header.txt MBS_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff MBS4527-1M-LZ.length.txt MBS4527-2M-LZ.length.txt
diff MBS4527-1M-LZ.length.txt MBS4527-3M-LZ.length.txt
sort MBS4527-1M-LZ.length.txt > MBS_LZ.length.txt

#MBS RS
ls ../MBS*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 MBS4527-1M-RS.tpm.txt) <(sort -k 1 MBS4527-2M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 MBS4527-3M-RS.tpm.txt) > MBS_RS.combined.tpm.txt
echo -e "\tMBS4527-1M-RS.tpm.txt\tMBS4527-2M-RS.tpm.txt\tMBS4527-3M-RS.tpm.txt" | cat - MBS_RS.combined.tpm.txt > MBS_RS.combined.tpm.header.txt
rm MBS_RS.combined.tpm.txt
mv MBS_RS.combined.tpm.header.txt MBS_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff MBS4527-1M-RS.length.txt MBS4527-2M-RS.length.txt
diff MBS4527-1M-RS.length.txt MBS4527-3M-RS.length.txt
sort MBS4527-1M-RS.length.txt > MBS_RS.length.txt

#PPPP LZ
ls ../PPPP*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 PPPP98-3M-LZ.tpm.txt) <(sort -k 1 PPPP98-4M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 PPPP98-5M-LZ.tpm.txt) > PPPP_LZ.combined.tpm.txt
echo -e "\tPPPP98-3M-LZ.tpm.txt\tPPPP98-4M-LZ.tpm.txt\tPPPP98-3M-LZ.tpm.txt" | cat - PPPP_LZ.combined.tpm.txt > PPPP_LZ.combined.tpm.header.txt
rm PPPP_LZ.combined.tpm.txt
mv PPPP_LZ.combined.tpm.header.txt PPPP_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff PPPP98-3M-LZ.length.txt PPPP98-4M-LZ.length.txt
diff PPPP98-3M-LZ.length.txt PPPP98-5M-LZ.length.txt
sort PPPP98-3M-LZ.length.txt > PPPP_LZ.length.txt

#PPPP RS
ls ../PPPP*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 PPPP98-3M-RS.tpm.txt) <(sort -k 1 PPPP98-4M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 PPPP98-5M-RS.tpm.txt) > PPPP_RS.combined.tpm.txt
echo -e "\tPPPP98-3M-RS.tpm.txt\tPPPP98-4M-RS.tpm.txt\tPPPP98-3M-RS.tpm.txt" | cat - PPPP_RS.combined.tpm.txt > PPPP_RS.combined.tpm.header.txt
rm PPPP_RS.combined.tpm.txt
mv PPPP_RS.combined.tpm.header.txt PPPP_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff PPPP98-3M-RS.length.txt PPPP98-4M-RS.length.txt
diff PPPP98-3M-RS.length.txt PPPP98-5M-RS.length.txt
sort PPPP98-3M-RS.length.txt > PPPP_RS.length.txt

#SEG LZ
ls ../SEG*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 SEG4130-LZ.tpm.txt) <(sort -k 1 SEG4156-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4176-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4197-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4700-1M-LZ.tpm.txt) > SEG_LZ.combined.tpm.txt
echo -e "\tSEG4130-LZ.tpm.txt\tSEG4156-LZ.tpm.txt\tSEG4176-LZ.tpm.txt\tSEG4197-LZ.tpm.txt\tSEG4700-1M-LZ.tpm.txt" | cat - SEG_LZ.combined.tpm.txt > SEG_LZ.combined.tpm.header.txt
rm SEG_LZ.combined.tpm.txt
mv SEG_LZ.combined.tpm.header.txt SEG_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff SEG4130-LZ.length.txt SEG4156-LZ.length.txt
diff SEG4130-LZ.length.txt SEG4176-LZ.length.txt
diff SEG4130-LZ.length.txt SEG4197-LZ.length.txt
diff SEG4130-LZ.length.txt SEG4700-1M-LZ.length.txt
sort SEG4130-LZ.length.txt > SEG_LZ.length.txt

#SEG RS
ls ../SEG*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 SEG4130-RS.tpm.txt) <(sort -k 1 SEG4156-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4176-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4197-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 SEG4700-1M-RS.tpm.txt) > SEG_RS.combined.tpm.txt
echo -e "\tSEG4130-RS.tpm.txt\tSEG4156-RS.tpm.txt\tSEG4176-RS.tpm.txt\tSEG4197-RS.tpm.txt\tSEG4700-1M-RS.tpm.txt" | cat - SEG_RS.combined.tpm.txt > SEG_RS.combined.tpm.header.txt
rm SEG_RS.combined.tpm.txt
mv SEG_RS.combined.tpm.header.txt SEG_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff SEG4130-RS.length.txt SEG4156-RS.length.txt
diff SEG4130-RS.length.txt SEG4176-RS.length.txt
diff SEG4130-RS.length.txt SEG4197-RS.length.txt
diff SEG4130-RS.length.txt SEG4700-1M-RS.length.txt
sort SEG4130-RS.length.txt > SEG_RS.length.txt

#SFM LZ
ls ../SFM*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 SFM4513-LZ.tpm.txt) <(sort -k 1 SFM4514-LZ.tpm.txt)> SFM_LZ.combined.tpm.txt
echo -e "\tSFM4513-LZ.tpm.txt\tSFM4514-LZ.tpm.txt" | cat - SFM_LZ.combined.tpm.txt > SFM_LZ.combined.tpm.header.txt
rm SFM_LZ.combined.tpm.txt
mv SFM_LZ.combined.tpm.header.txt SFM_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff SFM4513-LZ.length.txt SFM4514-LZ.length.txt
sort SFM4513-LZ.length.txt > SFM_LZ.length.txt

#SFM RS
ls ../SFM*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 SFM4513-RS.tpm.txt) <(sort -k 1 SFM4514-RS.tpm.txt)> SFM_RS.combined.tpm.txt
echo -e "\tSFM4513-RS.tpm.txt\tSFM4514-RS.tpm.txt" | cat - SFM_RS.combined.tpm.txt > SFM_RS.combined.tpm.header.txt
rm SFM_RS.combined.tpm.txt
mv SFM_RS.combined.tpm.header.txt SFM_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff SFM4513-RS.length.txt SFM4514-RS.length.txt
sort SFM4513-RS.length.txt > SFM_RS.length.txt

#STF LZ
ls ../STF*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 STF4495-1M-LZ.tpm.txt) <(sort -k 1 STF4515-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 STF4516-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 STF4517-LZ.tpm.txt) > STF_LZ.combined.tpm.txt
echo -e "\tSTF4495-1M-LZ.tpm.txt\tSTF4515-LZ.tpm.txt\tSTF4516-LZ.tpm.txt\tSTF4517-LZ.tpm.txt" | cat - STF_LZ.combined.tpm.txt > STF_LZ.combined.tpm.header.txt
rm STF_LZ.combined.tpm.txt
mv STF_LZ.combined.tpm.header.txt STF_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff STF4495-1M-LZ.length.txt STF4515-LZ.length.txt
diff STF4495-1M-LZ.length.txt STF4516-LZ.length.txt
diff STF4495-1M-LZ.length.txt STF4517-LZ.length.txt
sort STF4495-1M-LZ.length.txt > STF_LZ.length.txt

#STF RS
ls ../STF*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 STF4495-1M-RS.tpm.txt) <(sort -k 1 STF4515-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 STF4516-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 STF4517-RS.tpm.txt) > STF_RS.combined.tpm.txt
echo -e "\tSTF4495-1M-RS.tpm.txt\tSTF4515-RS.tpm.txt\tSTF4516-RS.tpm.txt\tSTF4517-RS.tpm.txt" | cat - STF_RS.combined.tpm.txt > STF_RS.combined.tpm.header.txt
rm STF_RS.combined.tpm.txt
mv STF_RS.combined.tpm.header.txt STF_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff STF4495-1M-RS.length.txt STF4515-RS.length.txt
diff STF4495-1M-RS.length.txt STF4516-RS.length.txt
diff STF4495-1M-RS.length.txt STF4517-RS.length.txt
sort STF4495-1M-RS.length.txt > STF_RS.length.txt

#PAHNew LZ
ls ../PAHNew*LZ_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 PAHNew-1M-LZ.tpm.txt) <(sort -k 1 PAHNew-2M-LZ.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 PAHNew-4M-LZ.tpm.txt) > PAHNew_LZ.combined.tpm.txt
echo -e "\tPAHNew-1M-LZ.tpm.txt\tPAHNew-2M-LZ.tpm.txt\tPAHNew-4M-LZ.tpm.txt" | cat - PAHNew_LZ.combined.tpm.txt > PAHNew_LZ.combined.tpm.header.txt
rm PAHNew_LZ.combined.tpm.txt
mv PAHNew_LZ.combined.tpm.header.txt PAHNew_LZ.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff PAHNew-1M-LZ.length.txt PAHNew-2M-LZ.length.txt
diff PAHNew-1M-LZ.length.txt PAHNew-4M-LZ.length.txt
sort PAHNew-1M-LZ.length.txt > PAHNew_LZ.length.txt

#PAHNew RS
ls ../PAHNew*RS_counts.orthoIDappended.txt | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "_" -f 1)
        echo "${name}"
        Rscript 01_calculate_tpm.r "${file}"
done
#combine tpm output files for same genotype samples
join -1 1 -2 1 <(sort -k 1 PAHNew-1M-RS.tpm.txt) <(sort -k 1 PAHNew-2M-RS.tpm.txt) | join -1 1 -2 1 - <(sort -k 1 PAHNew-4M-RS.tpm.txt) > PAHNew_RS.combined.tpm.txt
echo -e "\tPAHNew-1M-RS.tpm.txt\tPAHNew-2M-RS.tpm.txt\tPAHNew-4M-RS.tpm.txt" | cat - PAHNew_RS.combined.tpm.txt > PAHNew_RS.combined.tpm.header.txt
rm PAHNew_RS.combined.tpm.txt
mv PAHNew_RS.combined.tpm.header.txt PAHNew_RS.combined.tpm.txt
#check to make sure gene lengths the same across all samples (check calculate_tpm_output.txt for these results)
diff PAHNew-1M-RS.length.txt PAHNew-2M-RS.length.txt
diff PAHNew-1M-RS.length.txt PAHNew-4M-RS.length.txt
sort PAHNew-1M-RS.length.txt > PAHNew_RS.length.txt


echo "Done!"
