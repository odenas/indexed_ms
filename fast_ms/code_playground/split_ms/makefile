MS = rep_100000000s_dis_500000t_abcd.t.ms

all : ../bin/dump_ms_blocks.x $(MS)
	$< -ms_path $(MS)  -block_size  14400 -out_prefix ciao -out_suffix .ms

#$(MS) : ../bin/matching_stats_parallel.x m.s m.t
#	$< -load_cst 1 -s_path m.s -t_path m.t

#clean : 
#	rm -f $(MS) $(MS).0 $(MS).0.range $(MS).4.range $(MS).rrr $(RIDX)
#	rm -f m.t.runs m.t.runs.0
