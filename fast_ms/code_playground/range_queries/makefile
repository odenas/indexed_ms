BS = 4

all : m.t.ms.$(BS).ridx
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 0 -to_idx 1 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 0 -to_idx 2 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 0 -to_idx 3 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 0 -to_idx 4 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 0 -to_idx 10 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 1 -to_idx 2 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 1 -to_idx 3 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 1 -to_idx 4 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 1 -to_idx 10 -block_size $(BS) -check 1
	../../bin/range_queries.x -ms_path m.t.ms -ridx_path $< -from_idx 8 -to_idx 10 -block_size $(BS) -check 1

m.t.ms.$(BS).ridx : ../../bin/dump_range_index.x m.t.ms
	$< -ms_path m.t.ms  -block_size $(BS)

m.t.ms : ../../bin/matching_stats_parallel.x m.s m.t
	$< -s_path m.s -t_path m.t 

clean : 
	rm -f m.t.ms m.t.ms.0 m.t.ms.0.range m.t.ms.4.range
	rm -f m.t.runs m.t.runs.0
