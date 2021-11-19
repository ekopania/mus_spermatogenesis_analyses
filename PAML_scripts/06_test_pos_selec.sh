#!/bin/bash

#PURPOSE: use a LRT to compare lnL values of M8 and M8a, use the chi2 function of PAML to test for significant LRT result (pos selection)

#rm pos_selec.txt
ls -d PAML_RESULTS/CODEML_OUTPUT_ENSMUSG00000072663 | while read dir;
#ls -d PAML_RESULTS/CODEML_OUTPUT* | head | while read dir;
do
	gene=$(echo "${dir}" | cut -d "/" -f 2 | cut -d "_" -f 3)
	echo "${gene}"
	lnL=$(grep "lnL" "${dir}/mlc" | tail -1 | awk '{print $5}')
	echo "${lnL}"
	lnLa=$(grep "lnL" "PAML_RESULTS_M8a/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}')
	echo "${lnLa}"
#	LRT=$(awk 'function abs(v) {return v < 0 ? -v : v} {print abs(2*($1-$2))}' <<< "$lnLa $lnL")
#	LRT=$(awk 'function abs(v) {return v < 0 ? -v : v} {print 2*($1-$2)}' <<< "$lnLa $lnL")
	LRT=$(awk 'function abs(v) {return v < 0 ? -v : v} {print 2*($2-$1)}' <<< "$lnLa $lnL")
	echo "${LRT}"
	p=$(chi2 1 "${LRT}" | awk '{print $8}')
	echo "${p}"
#	printf "%s\t%f\t%f\t%f\t%f\n" "${gene}" "${lnL}" "${lnLa}" "${LRT}" "${p}" >> pos_selec.txt
#	printf "%s\t%f\t%f\t%f\t%f\n" "${gene}" "${lnL}" "${lnLa}" "${LRT}" "${p}" >> pos_selec.noAbs.txt
#	printf "%s\t%f\t%f\t%f\t%f\n" "${gene}" "${lnL}" "${lnLa}" "${LRT}" "${p}" >> pos_selec.noAbs.reversed.txt
	printf "%s\t%f\t%f\t%f\t%f\n" "${gene}" "${lnL}" "${lnLa}" "${LRT}" "${p}" >> pos_selec.ENSMUSG00000072663.txt
done
