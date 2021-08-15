#!/bin/bash


TABLES=../../bin/diff_none_tables.x
FROM_THRESHOLD=8  # Threshold on MS values, as above.
TO_THRESHOLD=28  # The program builds tables for all such values of the threshold
MAX_N=300   # Max nZeros and nOnes
NEGATIVE=0
SAVE_ALL=0  # If 1, stores the full DAG rather than just the tables
TABLES_DIR=./tables

${TABLES} ${FROM_THRESHOLD} ${TO_THRESHOLD} ${MAX_N} ${MAX_N} ${NEGATIVE} ${SAVE_ALL} ${TABLES_DIR}

