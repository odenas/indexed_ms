#!/bin/bash

set -e

for d in $*;
do
	if [ -d $d ]; then
		echo "Testing $d ..."
		cd $d/test
		snakemake --cores 1 --delete-all-output && snakemake -p --cores 1
		if [ $? -ne 0 ]; then
			exit 1
		fi
		cd ../../;
		echo -e "DONE ($d) \n\n\n"
	else
		echo "Argument '$d' is not a directory. Skipping ..."
	fi;
done
