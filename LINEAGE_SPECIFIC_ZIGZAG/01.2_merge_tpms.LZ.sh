#!/bin/bash
#PURPOSE: Merge tpm files for different strains of the same (sub)species
#
# Job name:
#SBATCH --job-name=merge_tpm_LZ
#SBATCH --output=merge_tpm_output.LZ.txt
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#DOM
stuff1=$(head -1 BIK_LZ.combined.tpm.txt)
stuff2=$(head -1 DGA_LZ.combined.tpm.txt)
stuff3=$(head -1 LLLL_LZ.combined.tpm.txt)
stuff4=$(head -1 WWWW_LZ.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 $stuff4 > DOM_LZ.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 BIK_LZ.combined.tpm.txt) <(tail -n+2 DGA_LZ.combined.tpm.txt) > DOM_LZ.temp1.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 LLLL_LZ.combined.tpm.txt) <(tail -n+2 WWWW_LZ.combined.tpm.txt) > DOM_LZ.temp2.combined.tpm.txt
join -1 1 -2 1 DOM_LZ.temp1.combined.tpm.txt DOM_LZ.temp2.combined.tpm.txt >> DOM_LZ.spec.combined.tpm.txt
rm DOM_LZ.temp?.combined.tpm.txt

#MUS
stuff1=$(head -1 CCCC_LZ.combined.tpm.txt)
stuff2=$(head -1 MBS_LZ.combined.tpm.txt)
stuff3=$(head -1 PPPP_LZ.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 > MUS_LZ.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 CCCC_LZ.combined.tpm.txt) <(tail -n+2 MBS_LZ.combined.tpm.txt) > MUS_LZ.temp1.combined.tpm.txt
join -1 1 -2 1 MUS_LZ.temp1.combined.tpm.txt <(tail -n+2 PPPP_LZ.combined.tpm.txt) >> MUS_LZ.spec.combined.tpm.txt
rm MUS_LZ.temp1.combined.tpm.txt

#SPR
stuff1=$(head -1 SEG_LZ.combined.tpm.txt)
stuff2=$(head -1 SFM_LZ.combined.tpm.txt)
stuff3=$(head -1 STF_LZ.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 > SPR_LZ.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 SEG_LZ.combined.tpm.txt) <(tail -n+2 SFM_LZ.combined.tpm.txt) > SPR_LZ.temp1.combined.tpm.txt
join -1 1 -2 1 SPR_LZ.temp1.combined.tpm.txt <(tail -n+2 STF_LZ.combined.tpm.txt) >> SPR_LZ.spec.combined.tpm.txt
rm SPR_LZ.temp1.combined.tpm.txt

echo "Done!"
