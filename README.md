Fast matching statistics in small space
=========

Tools for computing a compact representation of the matching statistics array (denoted by MS in what follows) between an indexed text and a query string on any byte alphabet. The matching statistics array takes 2*m* bits (rather than *m* integers), where *m* is the length of the query. The program takes *O*(*m log c*) time and *O*(*n log c*) bits of memory, where *n* is the length of the text and *c* is the size of its alphabet. The query string is read twice in opposite directions, and it must be known offline; however, it is streamed from disk and never kept fully in memory. The MS array is streamed to disk as well.

Please open an issue if you have any problem running the tools. For longer questions about the algorithm or the experiments, you are also welcome to [send an email to Olgert](mailto:gertidenas@gmail.com).


References
------------

The theory behind this code, as well as a partial experimental study, are described in the following papers:

* D. Belazzougui, F. Cunial, and O. Denas (2018). [Fast matching statistics in small space.](https://drops.dagstuhl.de/opus/volltexte/2018/8952/) 17th International Symposium on Experimental Algorithms (SEA 2018). Schloss Dagstuhl-Leibniz-Zentrum fuer Informatik. [[SLIDES]](https://www.slideshare.net/FabioCunial/fast-matching-statistics-in-small-space)
* D. Belazzougui and F. Cunial (2014). [Indexed matching statistics and shortest unique substrings.](https://link.springer.com/chapter/10.1007/978-3-319-11918-2_18) International Symposium on String Processing and Information Retrieval. Springer. [[SLIDES]](https://www.slideshare.net/FabioCunial/indexed-matching-statistics-and-shortest-unique-substrings)

Please cite the SEA paper if you use this tool.


Requirements
------------

* A modern, C++11 ready compiler such as [g++](https://gcc.gnu.org) version 4.9 or higher, or [clang](https://clang.llvm.org) version 3.2 or higher.
* The [cmake](http://www.cmake.org) build system.
* A 64-bit operating system. The code was tested on both macOS and Linux.
* For reproducing the experiments: a standard installation of [RStudio](https://www.rstudio.com) with [Tidyverse](https://www.tidyverse.org/), and Python 3.


Installing and testing
------------

Clone this repository with the `--recursive` flag and `cd` to `indexed_ms`.  
Next, install the `sdsl-lite` requirements and `sdsl-lite` itself, with the command:

```
cd sdsl-lite
sh install.sh `pwd`/build
cd ..
```

this will create libraries used by our executables in `./sdsl-lite/build`. To build the executables of this package, run:

```
~$ cd fast_ms
~$ make
~$ cd ..
```

this should build the programs in `./fast_ms/bin` (files ending in `.x`). Run a program without arguments to get a description
of what it does and how to use.

To make it easy to run experiements and organize work we support [snakemake](https://snakemake.readthedocs.io/en/stable/index.html). To
get started, create a python virtual environment with Python 3.7 and install dependencies. E.g.

```
~$ conda create --prefix ./myenv python=3.7
~$ conda activate ./myenv
~$ pip install -r requirements.txt
```

With the virtual environment set up and active, make the documentation and example workflows

```
~$ make -C fast_ms/docs/ html
```

to view, point a browser to `fast_ms/docs/_build/html/index.html`.

With the environment active, you can also run tests:

```
make -C fast_ms/tests && make -C fast_ms/wrappers
```


Organization
------------

 - `sdsl-lite`: fork of `sdsl-lite` with optimizations on several operations on the suffix tree.
 - `fast_ms`: source and executables linked to our fork of `sdsl-lite`
 - `py_interface`: various python scripts and utilities 
 - `experiments`: the experiments for the manuscript 
 
 
Usage
------------



Related code
---------

The following software computes some form of matching statistics as well:

* 

