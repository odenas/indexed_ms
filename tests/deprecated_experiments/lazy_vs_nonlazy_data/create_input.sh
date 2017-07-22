#!/bin/bash


DUMP_CST=/Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/dump_cst
LENS=1000000000
LENT=5000000
FTEMPLATE=1Gs_5Mt

echo "generate rnd_${FTEMPLATE} ..."
python ../generate_input.py random abcde input_data/rnd_$FTEMPLATE.s $LENS  input_data/rnd_$FTEMPLATE.t $LENT
echo "fwd stree ..."
$DUMP_CST input_data/rnd_${FTEMPLATE}.s 
echo "bwd stree ..."
$DUMP_CST input_data/rnd_${FTEMPLATE}.t

for MP in 10 100 1000
do
	echo "generate ${FTEMPLATE}_$MP ..."
	python ../generate_input.py mutation abcdefghijklmnopqrst input_data/mut_${FTEMPLATE}_$MP.s $LENS  input_data/mut_${FTEMPLATE}_$MP.t $LENT --mutation_period $MP
	echo "fwd stree ..."
	$DUMP_CST input_data/mut_${FTEMPLATE}_$MP.s 
	echo "bwd stree ..."
	$DUMP_CST input_data/mut_${FTEMPLATE}_$MP.t
done

