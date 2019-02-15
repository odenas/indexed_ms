#!/bin/bash

S=rep_1.s
T=rep_1.t

mstat_index_input.py --num_blocks 4 --num_mut 1 rep 12 aaaaaaaaaab $S
../input_stats_data/exploration/index_based_query/node_iterator.x -load_cst 1 -load_maxrep 1 -min_node_depth 2 -max_str_depth 200 -dump_size 5 -s_path $S >$T
mstat_fast_ms.py --label lca --lca_parents --repeat 1 --load_cst $S $T
