#!/bin/bash

#PURPOSE: extract omega (dN/dS) for all genes, LZ expressed genes, RS genes; also saves LZ and RS expressed genes that don't have dN/dS estimates (see GENES_NOT_INCLUDED)

#All commands append to text files, so make sure they aren't already in directory
rm omega_list_RS.txt
rm RS_expressed_no_paml.txt
rm omega_list_RSinduced.txt
rm RS_induced_no_paml.txt
rm omega_list_LZ.txt
rm LZ_expressed_no_paml.txt
rm omega_list_LZinduced.txt
rm LZ_induced_no_paml.txt
rm omega_list.txt

cat ../../mus_expression_analysis/MULTI_MAP/gene_list_eve_RS_edgeR_wholeGenome.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_expressed_no_paml.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RS.txt; done

cat ../../mus_expression_analysis/MULTI_MAP/gene_list_RSinduced_edgeR_wholeGenome.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> RS_induced_no_paml.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_RSinduced.txt; done

cat ../../mus_expression_analysis/MULTI_MAP/gene_list_eve_LZ_edgeR_wholeGenome.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_expressed_no_paml.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZ.txt; done

cat ../../mus_expression_analysis/MULTI_MAP/gene_list_LZinduced_edgeR_wholeGenome.txt | while read gene; do ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" 2>> LZ_induced_no_paml.txt | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list_LZinduced.txt; done

ls -d PAML_RESULTS/CODEML_OUTPUT_ENSMUSG00000* | while read dir; do gene=$(echo "${dir}" | cut -d "_" -f 4); ds=$(grep "tree length for dS" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); dn=$(grep "tree length for dN" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $5}'); grep "omega" "PAML_RESULTS/CODEML_OUTPUT_${gene}/mlc" | awk '{print $4}' | sed "s/^/${gene} ${dn} ${ds} /" >> omega_list.txt; done

