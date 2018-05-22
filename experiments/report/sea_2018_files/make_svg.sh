#!/bin/bash

for p in figure-latex/*pdf
do
	echo $p;
	pdf2svg $p svg-figures/`basename $p .pdf`.svg
done 
