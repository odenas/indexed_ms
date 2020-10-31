# Introduction
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

# Example
## Simple run
Compilation should have produced a bunch of `.x` files under the `bin` directory like so.
```
odenas@tachylite-d:fast_ms$ ls bin/
compress_ms.o  dump_edge_lists.o  dump_nwd_lists.o    freq_profile.o  matching_stats_parallel.o  print_int_ms.o  range_queries.o          sandbox_maxrep.o       sandbox_selectatdist.o  split_ms.o
compress_ms.x  dump_edge_lists.x  dump_nwd_lists.x    freq_profile.x  matching_stats_parallel.x  print_int_ms.x  range_queries_profile.o  sandbox_maxrep.x       sandbox_selectatdist.x  split_ms.x
dump_cst.o     dump_maxrep.o      dump_range_index.o  input_stats.o   matching_stats_slow.o      print_ms.o      range_queries_profile.x  sandbox_parentcalls.o  split_compress_ms.o
dump_cst.x     dump_maxrep.x      dump_range_index.x  input_stats.x   matching_stats_slow.x      print_ms.x      range_queries.x          sandbox_parentcalls.x  split_compress_ms.x
```

To generate a matching statistics we will need to run 
`bin/matching_stats_parallel.x`. This program expects raw text files (just the sequence).
(We have plans to add FASTA support soon.) So, assuming target and query files are `input.txt` and `p1.txt` respectively, one can run
```
odenas@tachylite-d:fast_ms$ ./bin/matching_stats_parallel.x -s_path input.txt -t_path p1.txt  -nthreads 4
[...]
odenas@tachylite-d:fast_ms$ ls
input.txt  p1.txt  p1.txt.ms  p1.txt.ms.0  p1.txt.ms.1  p1.txt.ms.2  p1.txt.ms.3  p1.txt.runs  p1.txt.runs.0  p1.txt.runs.0_3  p1.txt.runs.1  p1.txt.runs.2  p1.txt.runs.3
odenas@tachylite-d:fast_ms$
```

The MS vector in binary form is the file `p1.txt.ms`. To get the int MS values one can run
```
odenas@tachylite-d:fast_ms$ ./bin/print_int_ms.x -ms_path p1.txt.ms
3 4 3 2 3 2 3 3 4 3 2 1
malloc_count ### exiting, total: 1,069,512, peak: 1,069,440, current: 1,024
odenas@tachylite-d:fast_ms$
```

To get the binary representation of, say the first 3 MS values, one can run
```
odenas@tachylite-d:fast_ms$ ./bin/print_ms.x -ms_path p1.txt.ms -start 0 -len 3
000
odenas@tachylite-d:fast_ms$ ./bin/print_ms.x -ms_path p1.txt.ms  -start 0 -len 24  # since p1.txt has length 12
000100111001100101001111
malloc_count ### exiting, total: 1,069,768, peak: 1,069,504, current: 1,024
odenas@tachylite-d:fast_ms$
```

One can also run max/sum range queries on the binary MS like so:
```
odenas@tachylite-d:fast_ms$ ./bin/range_queries.x -ms_path p1.txt.ms -from_idx 0 -to_idx 3 -compression none -algo t -op max
[0, 3): 4
odenas@tachylite-d:fast_ms$ ./bin/range_queries.x -ms_path p1.txt.ms -from_idx 0 -to_idx 3 -compression none -algo t -op sum
[0, 3): 10
odenas@tachylite-d:fast_ms$
```


## Repeated run
For large input sizes and/or repeated matching statistics computations on a target file, it saves a lot of time to compute partial results. The most important is the CST. In this case one can run
```
odenas@tachylite-d:fast_ms$  ./bin/dump_cst.x -s_path input.txt
 * loadding index string  from input.txt
 * building the CST of length 128 DONE (0seconds)
 * dumping the CST to input.txt.fwd.stree DONE (0seconds)
 * loadding index string  from input.txt
 * building the CST of length 128 DONE (0seconds)
 * dumping the CST to input.txt.rev.stree DONE (0seconds)
malloc_count ### exiting, total: 33,364,660, peak: 8,010,610, current: 0
odenas@tachylite-d:fast_ms$ ls
input.txt            input.txt.rev.stree  p1.txt.ms    p1.txt.ms.1  p1.txt.ms.3  p1.txt.runs.0    p1.txt.runs.1  p1.txt.runs.3
input.txt.fwd.stree  p1.txt               p1.txt.ms.0  p1.txt.ms.2  p1.txt.runs  p1.txt.runs.0_3  p1.txt.runs.2
```

and 2 `.stree` files were generated. To compute the MS now we run the same command + an extra command line option. (this should run faster).
```
odenas@tachylite-d:fast_ms$ ./bin/matching_stats_parallel.x -s_path input.txt -t_path p1.txt  -nthreads 4 -load_cst 1
```

The 2 files I used for the examples are
```
odenas@tachylite-d:fast_ms$ cat input.txt
bebbeabcaecdbddddaccdeacdebeeeaedbebdacecbdceaebbbbabcbacebabecdcebbaabdddaebbedbaccbdbadcaccdcdbeadabcbebddcdacdbeedcabbdcacdae
odenas@tachylite-d:fast_ms$ cat p1.txt
abbabaaabbaa
odenas@tachylite-d:fast_ms$
```
