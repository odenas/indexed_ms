#!/bin/bash


#for LENS in 100000 200000 400000 800000

#for LENS in 1000000
#do
#	python ../generate_input.py file ../proteins.50MB --base_dir input_data --len_s $LENS --len_t 100000
#done

for LENS in 10000000
do
	python ../generate_input.py random abcde --base_dir input_data --len_s $LENS --len_t 100000
done

#for MP in 10 20 40 80 160 320 640 1280 2560 5120 10240 20480 40960
#do
#	python ../generate_input.py mutation abcdefghijklmnopqrst --base_dir input_data --len_s 1000000 --len_t 100000 --mutation_period $MP
#done

