#!/bin/bash


#for LENS in 100000 200000 400000 800000

for LENS in 10000000 20000000 40000000 80000000
do
	python ../generate_input.py mutation acgt --base_dir input_data --mutation_period 1000 --len_s $LENS --len_t 200000
done


