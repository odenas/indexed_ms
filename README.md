<p align="left">
<img src="./logo.png" width="245" height="42"/><br>
<b>Fast and compact matching statistics</b>
</p>

Tools for computing, compressing, and analyzing a compact representation of the matching statistics (MS) array between an indexed text and a query string on any byte alphabet. Our representation of the MS array takes 2 bits per character (instead of one integer), and our compression programs can shrink it down to 0.8 or 0.2 bits per character (and in some cases to 0.04 bits per character), depending on the dataset.

Our sequential program takes *O*(*m log c*) time and *O*(*n log c*) bits of memory, where *n* is the length of the text, *m* is the length of the query, and *c* is the size of the alphabet. The query string is read twice in opposite directions, and it must be given offline; however, it is streamed from disk and never kept fully in memory. The MS array is streamed to disk as well. Several heuristics are available to speed up computation in practice (for example when the query is very similar to the text). Our parallel, shared-memory program is 30 times faster than our sequential program when using 48 cores.

We also provide fast tools for computing the maximum and the average MS value inside a range of positions.

Please open an issue if you have any problem running the tools. For longer questions about the algorithm or the experiments, you are also welcome to [send an email to Olgert](mailto:gertidenas@gmail.com).


References
------------

The theory behind this code, as well as a partial experimental study, are described in the following papers:

* F. Cunial, O. Denas, and D. Belazzougui (2021). Fast and compact matching statistics analytics. bioRxiv ________.
* D. Belazzougui, F. Cunial, and O. Denas (2018). [Fast matching statistics in small space.](https://drops.dagstuhl.de/opus/volltexte/2018/8952/) 17th International Symposium on Experimental Algorithms (SEA 2018). Schloss Dagstuhl-Leibniz-Zentrum fuer Informatik. [[SLIDES]](https://www.slideshare.net/FabioCunial/fast-matching-statistics-in-small-space)
* D. Belazzougui and F. Cunial (2014). [Indexed matching statistics and shortest unique substrings.](https://link.springer.com/chapter/10.1007/978-3-319-11918-2_18) International Symposium on String Processing and Information Retrieval. Springer. [[SLIDES]](https://www.slideshare.net/FabioCunial/indexed-matching-statistics-and-shortest-unique-substrings)

Please cite the SEA paper if you use this tool.


Requirements
------------

* A compiler that supports C++11 and OpenMP, such as [GCC](https://gcc.gnu.org) (we compiled successfully on versions 6.2, 8.3, 10.3).
* The [cmake](http://www.cmake.org) build system.
* A 64-bit operating system (we tested successfully on macOS 10.15, CentOS, Ubuntu).
* For reproducing the experiments: a standard installation of [RStudio](https://www.rstudio.com) with [Tidyverse](https://www.tidyverse.org/), and Python 3.


Installing and testing
------------

Clone this repository and `cd` to `indexed_ms`. Next, install the `sdsl-lite` requirements and `sdsl-lite` itself, with the command:

```
cd sdsl-lite; sh install.sh `pwd`/build ; cd ..
```

This will create libraries used by our executables in `./sdsl-lite/build`. To build the executables of this package, run:

```
cd fast_ms; make; cd ..
```

This should build the executables in `./fast_ms/bin` (files ending in `.x`). Run a program without arguments to get a description of what it does and of how to use it.

To make it easy to run experiements and organize work we support [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). To get started, create a Python virtual environment with Python 3.7 and install dependencies. E.g.

```
conda create --prefix ./myenv python=3.7
conda activate ./myenv
pip install -r requirements.txt
```

With the virtual environment set up and active, make the documentation and example workflows:

```
make -C fast_ms/docs/ html
```

To view the documentation, point a browser to `fast_ms/docs/_build/html/index.html`.

With the environment active, you can also run tests:

```
make -C fast_ms/tests && make -C fast_ms/wrappers
```


Directory structure
------------

 - `sdsl-lite`: fork of `sdsl-lite` with optimizations on several operations on the compressed suffix tree.
 - `fast_ms`: our main source code and executables based on our fork of `sdsl-lite`.
 - `py_interface`: various Python scripts and utilities.
 - `experiments`: the experiments for the manuscripts.


Usage
------------

After compilation, the `fast_ms/bin` directory should contain the following executable `.x` files:
```
$ ls ./fast_ms/bin/*.x
bin/compress_ms.x            bin/matching_stats_parallel.x
bin/diff_compress_ms.x       bin/matching_stats_slow.x
bin/diff_long_compress_ms.x  bin/matching_stats.x
bin/diff_long_tables.x       bin/print_freq.x
bin/diff_none_compress_ms.x  bin/print_int_ms.x
bin/diff_none_tables.x       bin/print_ms.x
bin/diff_tables.x            bin/range_queries_profile.x
bin/dump_cst.x               bin/range_queries.x
bin/dump_edge_lists.x        bin/sandbox_maxrep.x
bin/dump_maxrep.x            bin/sandbox_parentcalls.x
bin/dump_nwd_lists.x         bin/sandbox_selectatdist.x
bin/dump_range_index.x       bin/split_compress_ms.x
bin/freq_profile.x           bin/split_ms.x
bin/input_stats.x
```
Every one of these programs will display a help message if it is ran without parameters.

In this tutorial we will show how to build the matching statistics of a query `q1.txt` with respect to an index file `index.txt`. Then we will show how to execute range-max and range-sum queries on a matching statistics array.


The index and query files I will for the examples are
```
(myenv) user@laptop:fast_ms$ cat examples/tutorial_files/input.txt
bebbeabcaecdbddddaccdeacdebeeeaedbebdacecbdceaebbbbabcbacebabecdcebbaabdddaebbedbaccbdbadcaccdcdbeadabcbebddcdacdbeedcabbdcacdae
(myenv) user@laptop:fast_ms$ cat examples/tutorial_files/q1.txt
abbabaaabbaa
```
We will copy them to this directory for convenience
```
(myenv) user@laptop:fast_ms$ cp examples/tutorial_files/*.txt .
(myenv) user@laptop:fast_ms$ ls *txt
index.txt  q1.txt
```

The matching statistics algorithm relies on the suffix tree of the forward and backward index file. Since this file
is fixed, it makes sense to build this once and store it to disk.

```
(myenv) user@laptop:fast_ms$ bin/dump_cst.x -s_path index.txt
 * loadding index string  from index.txt 
 * building the CST of length 129 DONE (0seconds)
 * dumping the CST to index.txt.fwd.stree DONE (0seconds)
 * loadding index string  from index.txt 
 * building the CST of length 129 DONE (0seconds)
 * dumping the CST to index.txt.rev.stree DONE (0seconds)
malloc_count ### exiting, total: 33,375,722, peak: 8,010,898, current: 0
(myenv) user@laptop:fast_ms$ ls *stree
index.txt.fwd.stree  index.txt.rev.stree
```
another optional partial index that can speed up matching statistics is the maxrep vector containing
a list of nodes that maximal repeats.

```
(myenv) user@laptop:fast_ms$ bin/dump_maxrep.x -s_path index.txt -load_cst 1
 * loading the CST from index.txt.rev.stree DONE (0 seconds)
 * computing MAXREP DONE (0 milliseconds)
 * dumping MAXREP to index.txt.rev.maxrep DONE (0 seconds)
malloc_count ### exiting, total: 22,407, peak: 12,444, current: 0
(myenv) user@laptop:fast_ms$ ls *maxrep
index.txt.rev.maxrep
```

The flag `-load_cst 1` tells the program to load the CST file generated above. The name of the output file
suggests that we used the suffix tree of the reverse index.

Now we are ready to build the matching statistics type, for example,

```
(myenv) user@laptop:fast_ms$ bin/matching_stats_parallel.x -s_path index.txt -t_path q1.txt -load_cst 1 -load_maxrep 1 -lazy_wl 1 -nthreads 4
building RUNS ... 
 * loading the CST from index.txt.fwd.stree DONE (0 seconds, 130 leaves)
 ** filling 4 slices with : 4 threads ...
 *** [0 .. 3)
 *** [3 .. 6)
 *** [6 .. 9)
 *** [9 .. 12)
 ** DONE: 0
 ** correcting 1 intervals over 4 threads ... 
 *** [0][[10,1), node(24, 24, 45, 256, 46)]     intervals 0 - 3
 ** DONE: 0
 * merging into index.txt_q1.txt.runs ... 
DONE (0 seconds)
building MS ... 
 * loading the CST from index.txt.rev.stree DONE (0 seconds, 130 leaves)
 ** launching ms computation over : [0 .. 3)
 ** launching ms computation over : [3 .. 6)
 ** launching ms computation over : [6 .. 9)
 ** launching ms computation over : [9 .. 12)
 *** [0]filled 9 of 24 entries 
 *** [1]filled 7 of 24 entries 
 *** [2]filled 12 of 24 entries 
 *** [3]filled 6 of 24 entries 
 * merging into index.txt_q1.txt.ms ... 
 ** adding [0 .. 3) from index.txt_q1.txt.ms.0
 ** adding [3 .. 6) from index.txt_q1.txt.ms.1
 ** adding [6 .. 9) from index.txt_q1.txt.ms.2
 ** adding [9 .. 12) from index.txt_q1.txt.ms.3
DONE (0 seconds)
malloc_count ### exiting, total: 623,703, peak: 125,656, current: 1,571
(myenv) user@laptop:fast_ms$ ls *ms
index.txt_q1.txt.ms
```

The flags `-load_cst 1` and `load_maxrep 1` say that we want to load the suffix tree and maxrep files from disk instead of building them on the fly.
They are optional. You can also see that we used 4 threads to build the vector.
The matching statistics file name is of the form `<index file name>_<query file name>.ms`. The file 
is a binary bit vector. To get the int MS values one can run

```
(myenv) user@laptop:fast_ms$ bin/print_int_ms.x -ms_path index.txt_q1.txt.ms
3 4 3 2 3 2 3 3 4 3 2 1
malloc_count ### exiting, total: 1,069,512, peak: 1,069,440, current: 1,024
```

To get the binary representation of, say the first 3 MS values, one can run
```
(myenv) user@laptop:fast_ms$ bin/print_ms.x -ms_path index.txt_q1.txt.ms -start 0 -len 3
000
(myenv) user@laptop:fast_ms$ bin/print_ms.x -ms_path index.txt_q1.txt.ms -start 0 -len 24  # since q1.txt has length 12
000100111001100101001111
malloc_count ### exiting, total: 1,069,768, peak: 1,069,504, current: 1,024
```

We are ready to run range queries now. For example:

```
(myenv) user@laptop:fast_ms$ bin/range_queries.x -ms_path index.txt_q1.txt.ms -from_idx 1 -to_idx 2 -compression none -algo t -op max
[1, 2) 1: 4
```

We support range queries on compressed MS vectors, so in this case `-compression none` indicates that `-ms_path` is not compressed. The
`t` algorithm is slower than the `d`, but supports more formats. For a large number of queries it makes sense to compute partial
results.

```
(myenv) user@laptop:fast_ms$ bin/dump_range_index.x -ms_path index.txt_q1.txt.ms -block_size 2 -op max
(myenv) user@laptop:fast_ms$ ls *ridx
index.txt_q1.txt.ms.2.ridx
```
The file name is in the format `<ms file name>_<block size>.ridx`.
The block size depends on the size of your index. A small block size will take longer to compute
and result in a bigger index, but in turn it will provide a bigger speed up. Having blocks
sizes up to 10s of millions should yield good results.
Then for the queries, for example

```
(myenv) user@laptop:fast_ms$ bin/range_queries.x -ms_path index.txt_q1.txt.ms -ridx_path index.txt_q1.txt.ms.2.ridx -block_size 2 -from_idx 1 -to_idx 2 -compression none -algo t -op max
[1, 2) 1: 4
```


One can also run max/sum range queries on the binary MS like so:
```
(myenv) user@laptop:fast_ms$ ./bin/range_queries.x -ms_path q1.txt.ms -from_idx 0 -to_idx 3 -compression none -algo t -op max
[0, 3): 4
(myenv) user@laptop:fast_ms$ ./bin/range_queries.x -ms_path q1.txt.ms -from_idx 0 -to_idx 3 -compression none -algo t -op sum
[0, 3): 10
(myenv) user@laptop:fast_ms$
```

We also provide snakemake wrappers for the main utilities.
Documentation, if built as described in the tutorial in the parent
folder, is provided in `docs/_build/html/index.html`. A working
example that replicates this tutorial is provided in `examples/workflows`.





