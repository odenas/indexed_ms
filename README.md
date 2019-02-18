# Fast matching statistics in small space

Implementation and analysis of algorithms for the matching statistics problem. 

## Organization

 - `fast_ms`: source and executables linked to our fork of `sdsl-lite`
 - `py_interface`: python scripts used as an interface to executables in `fast_ms` and other utility scripts
 - `sdsl-lite`: fork of `sdsl-lite` with optimizations on several operations on the suffix tree.
 - `tests`: the experiments. 
 
## Usage
### Requirements

 - Requirements of `sdsl-lite`: https://github.com/simongog/sdsl-lite#requirements
 - Python version 3
 - (for reports only) A standard https://www.rstudio.com/ installation with [tidyverse](https://www.tidyverse.org/)

### Installation
Clone this repository with the `--recursive` flag and `cd` to `indexed_ms`.  
Next, install `sdsl-lite` requirements, then `sdsl-lite` itself with the command:

```
~$ cd sdsl-lite
~$ sh install.sh `pwd`/build
~$ cd ..
```

this will create libraries used by our executables in `./sdsl-lite/build`.
To build out executables run

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
