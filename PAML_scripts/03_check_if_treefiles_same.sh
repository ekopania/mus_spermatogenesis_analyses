#PURPOSE: I manually parallelized alignments and treefile generation for each gene and ended up with multiple CODEML_ directories for some genes; check to make sure these generated the same tree before deleting one version

ls -d ../TEMP*/CODEML_OUTPUT_ENSMUSG00000* | wc
dif_count=0
ls -d ../TEMP*/CODEML_OUTPUT_ENSMUSG00000* | while read dir
do
	name=$(echo "${dir}" | cut -d "/" -f 3)
	echo "${name}"
	num_diffs=$(diff "${dir}/final_alignment.phy.treefile" "${name}/final_alignment.phy.treefile" | wc -l)
	echo "${num_diffs}"
	if [ "${num_diffs}" -gt 0 ]; then
		dif_count=$((dif_count + 1))
	elif [ "${num_diffs}" -eq 0 ]; then
		rm -r "${dir}"
	fi
	echo "${dif_count}"
done
#echo "${dif_count}"
