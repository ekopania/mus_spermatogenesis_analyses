#!/bin/bash

#Purpose: Loop through coding genes to run codeml separately for each gene

ls -d CODEML_OUTPUT* | while read dir
do
        echo "${dir}"; cd "${dir}"; codeml ../codeml.ctl &> codeml.out < /dev/null; cd /scratch/1/kopania/PAML_runs/WHOLE_GENOMES/; echo "${dir}" &
done

