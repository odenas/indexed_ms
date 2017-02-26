#!/bin/bash


LENS=100000000
LENT=5000000

python ../generate_input.py random abcde input_data/rnd_10Ms_5Mt.s $LENS  input_data/rnd_1Ms_500kt.t $LENT

for MP in 10 100 1000
do
	python ../generate_input.py mutation abcdefghijklmnopqrst input_data/mut_10Ms_5Mt_$MP.s $LENS  input_data/mut_10Ms_5Mt_$MP.t $LENT --mutation_period $MP
done

