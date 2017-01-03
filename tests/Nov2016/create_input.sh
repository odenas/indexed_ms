#!/bin/bash


for LENS in 100000 150000 200000 250000
do
	for LENT in 51200
	do
		python ../generate_input.py non_random  ../proteins.50MB --base_dir input_data --len_s $LENS --len_t $LENT
	done
done


