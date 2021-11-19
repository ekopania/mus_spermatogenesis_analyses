#!/bin/bash
#PURPOSE: Merge tpm files for different strains of the same (sub)species
#
# Job name:
#SBATCH --job-name=merge_tpm_RS
#SBATCH --output=merge_tpm_output.RS.txt
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
stuff1=$(head -1 BIK_RS.combined.tpm.txt)
stuff2=$(head -1 DGA_RS.combined.tpm.txt)
stuff3=$(head -1 LLLL_RS.combined.tpm.txt)
stuff4=$(head -1 WWWW_RS.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 $stuff4 > DOM_RS.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 BIK_RS.combined.tpm.txt) <(tail -n+2 DGA_RS.combined.tpm.txt) > DOM_RS.temp1.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 LLLL_RS.combined.tpm.txt) <(tail -n+2 WWWW_RS.combined.tpm.txt) > DOM_RS.temp2.combined.tpm.txt
join -1 1 -2 1 DOM_RS.temp1.combined.tpm.txt DOM_RS.temp2.combined.tpm.txt >> DOM_RS.spec.combined.tpm.txt
rm DOM_RS.temp?.combined.tpm.txt

#MUS
stuff1=$(head -1 CCCC_RS.combined.tpm.txt)
stuff2=$(head -1 MBS_RS.combined.tpm.txt)
stuff3=$(head -1 PPPP_RS.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 > MUS_RS.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 CCCC_RS.combined.tpm.txt) <(tail -n+2 MBS_RS.combined.tpm.txt) > MUS_RS.temp1.combined.tpm.txt
join -1 1 -2 1 MUS_RS.temp1.combined.tpm.txt <(tail -n+2 PPPP_RS.combined.tpm.txt) >> MUS_RS.spec.combined.tpm.txt
rm MUS_RS.temp1.combined.tpm.txt

#SPR
stuff1=$(head -1 SEG_RS.combined.tpm.txt)
stuff2=$(head -1 SFM_RS.combined.tpm.txt)
stuff3=$(head -1 STF_RS.combined.tpm.txt)
echo $stuff1 $stuff2 $stuff3 > SPR_RS.spec.combined.tpm.txt
join -1 1 -2 1 <(tail -n+2 SEG_RS.combined.tpm.txt) <(tail -n+2 SFM_RS.combined.tpm.txt) > SPR_RS.temp1.combined.tpm.txt
join -1 1 -2 1 SPR_RS.temp1.combined.tpm.txt <(tail -n+2 STF_RS.combined.tpm.txt) >> SPR_RS.spec.combined.tpm.txt
rm SPR_RS.temp1.combined.tpm.txt

echo "Done!"
