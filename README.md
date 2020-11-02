# Fast matching statistics in small space

Tools for computing a compact representation of the matching statistics array (denoted by MS in what follows) between an indexed text and a query string on any byte alphabet. The matching statistics array takes 2*m* bits (rather than *m* integers), where *m* is the length of the query. The program takes *O*(*m log c*) time and *O*(*n log c*) bits of memory, where *n* is the length of the text and *c* is the size of its alphabet. The query string is read twice in opposite directions, and it must be known offline; however, it is streamed from disk and never kept fully in memory. The MS array is streamed to disk as well.

## References

The theory behind this code, as well as a partial experimental study, are described in the following papers:

* D. Belazzougui, F. Cunial, and O. Denas (2018). [Fast matching statistics in small space.](https://drops.dagstuhl.de/opus/volltexte/2018/8952/) 17th International Symposium on Experimental Algorithms (SEA 2018). Schloss Dagstuhl-Leibniz-Zentrum fuer Informatik. [[SLIDES]](https://www.slideshare.net/FabioCunial/fast-matching-statistics-in-small-space)
* D. Belazzougui and F. Cunial (2014). [Indexed matching statistics and shortest unique substrings.](https://link.springer.com/chapter/10.1007/978-3-319-11918-2_18) International Symposium on String Processing and Information Retrieval. Springer. [[SLIDES]](https://www.slideshare.net/FabioCunial/indexed-matching-statistics-and-shortest-unique-substrings)

Please cite the SEA paper if you use this tool.

## Organization

 - `sdsl-lite`: fork of `sdsl-lite` with optimizations on several operations on the suffix tree.
 - `fast_ms`: source and executables linked to our fork of `sdsl-lite`
 - `py_interface`: various python scripts and utilities 
 - `experiments`: the experiments for the manuscript 
 
## Usage
### Requirements

 - Requirements of `sdsl-lite`: https://github.com/simongog/sdsl-lite#requirements
 - (for experiments only) A standard https://www.rstudio.com/ installation with [tidyverse](https://www.tidyverse.org/) and Python version 3

### Installation
Clone this repository with the `--recursive` flag and `cd` to `indexed_ms`.  
Next, install `sdsl-lite` requirements, then `sdsl-lite` itself with the command:

```
~$ cd sdsl-lite
~$ sh install.sh `pwd`/build
~$ cd ..
```

this will create libraries used by our executables in `./sdsl-lite/build`. To build out executables run

```
~$ cd fast_ms
~$ make
~$ cd tests
~$ make  # this will run tests  
~$ cd ../../
```

this should build the programs in `./fast_ms/bin`. Run a program without arguments to get a description
of what it does and how to use.


# Contact
For suggestions or questions please write to Olgert Denas `gertidenas@gmail.com`
