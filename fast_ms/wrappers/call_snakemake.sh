#!/bin/bash


for d in $*;
do
	if [ -d $d ]; then
		echo "Testing $d ..."
		cd $d/test
		snakemake --cores 1 --delete-all-output && snakemake -p --cores 1
		cd ../../;
		echo -e "DONE ($d) \n\n\n"
	else
		echo "Argument '$d' is not a directory. Skipping ..."
	fi;
done
