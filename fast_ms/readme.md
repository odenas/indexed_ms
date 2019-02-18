This directory contains the `fast_ms` executables. Follow the instructions in 
the parent directory, then type `make` to build them. 

A typical workflow starts with an index file (`index.txt`), and one or more queries (`q1.txt, ..., qk.txt`). 
Since, most likely, the programs below will be called multiple times, it makes sense to build the
indexes first. 

 - The CST: `bin/dump_cst -s_path index.txt` will create the CST for the forward and backward index.
 - The maxrep vector containing a list of nodes that are maximal repeats. `bin/dump_maxrep.x -s_path index.txt -load_cst 1` 
   will create the binary file.

To build the matching statistics type, for example, 
`bin/matching_stats_parallel.x -s_path index.txt -t_path qi.txt -load_cst 1 -load_maxrep 1 -lazy_wl 1 -nthreads 4` 

To run range queries, first build a range query index. For example:
`bin/dump_range_index.x -ms_path index.txt.ms -block_size 1024`
Then for the queries, for example 
`bin/range_queries.x -ms_path index.txt.ms -ridx_path index.txt.ms.1024.ridx -block_size 1024 -from_idx 1 -to_idx 5`

Other executables are used for various experimental features and analysis. 
