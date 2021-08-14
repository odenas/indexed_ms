#!/bin/bash


TABLES=../../../bin/diff_none_tables.x
COMPRESS=../../../bin/diff_none_compress_ms.x

THRESHOLD=5  # Threshold on MS values, as above.
NEGATIVE=1
GREEDY=0  # Greedy strategy for too large windows. Either 0 or 1.
N_ZEROS=300  # Max n. of zeros in a precomputed table
N_ONES=300  # Max n. of ones in a precomputed table

TABLES_DIR=./tables
OUTPUT_DIR=odir
LOG_FILE=log_file.txt
MS_FILE=inp.s_inp.t.ms

${COMPRESS} -ms_path ${MS_FILE} -compression rle -outputDir ${OUTPUT_DIR} -tablesDir ${TABLES_DIR} \
	-threshold ${THRESHOLD} -nZeros ${N_ZEROS} -nOnes ${N_ONES} -negative ${NEGATIVE} -greedy ${GREEDY} -verbose 0  1> ${LOG_FILE} 2>&1

