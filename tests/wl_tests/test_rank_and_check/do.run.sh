#!/bin/bash



#for fname in rank_and_check__uv single_rank__uv rank_and_check__has_wl single_rank__has_wl
for fname in ../../xcode_bin/fd_ms
do
	echo "*************************************"
	date
	echo "*************************************"
	EP=$fname make -n
	make clean 
	echo "*************************************"
	date
	echo "*************************************"
done

