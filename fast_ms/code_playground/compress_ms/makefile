BS = 4
MS = m.t.ms
RIDX = $(MS).$(BS).ridx

all : ../../bin/range_queries.x $(RIDX) $(MS).rrr
	$< -ms_path $(MS).rrr -ridx_path $(RIDX) -from_idx 0 -to_idx 1 -block_size $(BS) -check 1

$(MS).rrr : ../../bin/compress_ms.x $(MS)
	$< -ms_path $(MS)

$(RIDX) : ../../bin/dump_range_index.x $(MS)
	$< -ms_path $(MS)  -block_size $(BS)

$(MS) : ../../bin/matching_stats_parallel.x m.s m.t
	$< -load_cst 1 -s_path m.s -t_path m.t

clean : 
	rm -f $(MS) $(MS).0 $(MS).0.range $(MS).4.range $(MS).rrr $(RIDX)
	rm -f m.t.runs m.t.runs.0
