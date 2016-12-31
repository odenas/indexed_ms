#!/bin/bash


for LENS in 400 800 1600 3200 6400 12800 25600
do
	for LENT in 400 800 1600 3200 6400 12800 25600
	do
		python ../generate_input.py non_random  ../proteins.50MB --base_dir input_data --len_s $LENS --len_t $LENT
	done
done


