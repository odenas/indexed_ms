BS = 5
IID = mut_200s_64t_15

RQ_ARG = -compression none -block_size $(BS)

all : $(IID).t.ms.$(BS).ridx
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 0 -to_idx 1 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 0 -to_idx 2 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 0 -to_idx 3 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 0 -to_idx 4 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 0 -to_idx 10 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 1 -to_idx 2 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 1 -to_idx 3 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 1 -to_idx 4 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 1 -to_idx 10 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 8 -to_idx 10 -check 1
	../../bin/range_queries.x $(RQ_ARG) -ms_path $(IID).t.ms -ridx_path $< -from_idx 25 -to_idx 26 -check 1

$(IID).t.ms.$(BS).ridx : ../../bin/dump_range_index.x $(IID).t.ms
	$< -op max -ms_path $(IID).t.ms  -block_size $(BS)


$(IID).t.ms : ../../bin/matching_stats_parallel.x $(IID).s $(IID).t
	$< -s_path $(IID).s -t_path $(IID).t 

clean : 
	rm -f $(IID).t.ms $(IID).t.ms.0 $(IID).t.ms.0.range $(IID).t.ms.4.range
	rm -f $(IID).t.ms.$(BS).ridx
	rm -f $(IID).t.runs $(IID).t.runs.0
