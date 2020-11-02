# Fast matching statistics in small space

Tools for computing a compact representation of the matching statistics array (denoted by MS in what follows) between an indexed text and a query string on any byte alphabet. The matching statistics array takes 2*m* bits (rather than *m* integers), where *m* is the length of the query. The program takes *O*(*n log c*) bits of memory, where *n* is the length of the text and *c* is the size of its alphabet, and *O*(*m log c*) time. It consists in two passes over the query string, which must be given offline but is streamed from disk and never kept fully in memory.

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
