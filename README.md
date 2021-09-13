<p>
<img align="center" src="./logo.png" width="245" height="42"/><br>
<b>Fast and compact matching statistics</b>
</p>

Toolkit for computing and analyzing a compact representation of the matching statistics array (denoted by MS in what follows) between an indexed text and a query string on any byte alphabet. Our encoding of the matching statistics array takes 2*m* bits (rather than *m* integers), where *m* is the length of the query. Our sequential program takes *O*(*m log c*) time and *O*(*n log c*) bits of memory, where *n* is the length of the text and *c* is the size of its alphabet. The query string is read twice in opposite directions, and it must be given offline; however, it is streamed from disk and never kept fully in memory. The MS array is streamed to disk as well.

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
* A 64-bit operating system (we tested the code on macOS 10.15, CentOS, Ubuntu).
* For reproducing the experiments: a standard installation of [RStudio](https://www.rstudio.com) with [Tidyverse](https://www.tidyverse.org/), and Python 3.


Installing and testing
------------

Clone this repository with the `--recursive` flag and `cd` to `indexed_ms`.  
Next, install the `sdsl-lite` requirements and `sdsl-lite` itself, with the command:

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
