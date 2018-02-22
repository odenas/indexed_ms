# Experiments for "Fast matching statistics in small space" 


## Organization

 - `fast_ms`: source and executables linked to our fork of `sdsl-lite`
 - `py_interface`: python scripts used as an interface to executables in `fast_ms` and other utility scripts
 - `sdsl-lite`: fork of `sdsl-lite` with optimizations on several operations on the suffix tree.
 - `tests`: the experiments. 
 
## Usage
### Requirements

 - Requirements of `sdsl-lite`: https://github.com/simongog/sdsl-lite#requirements
 - Python version 3
 - (for reports only) A standard https://www.rstudio.com/ installation with `tidyverse`

### Installation
This installation guide has been tested on OSX and linux. It should be possible to use in windows as well, but we haven't tested it. First, clone this repository with the `--recursive` flag (this tells git to pull all submodules) in a working directory, say `~/matching_statistics`. 

Next, install `sdsl-lite` by running the following commands:

```
~$ cd ~/matching_statistics/indexed_ms/sdsl-lite
~$ bash install.sh `pwd`/build
```

this will create libraries used by our executables in `~/matching_statistics/indexed_ms/sdsl-lite/build`. Now we build our executables in `fast_ms` as follows. Edit `~/matching_statistics/indexed_ms/fast_ms/makefile` changing the lines 

```
SDSL_BASE_DIR = /home/brt/code/matching_statistics/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include    
INCLUDES_PRIME = /home/cunial/arch/Darwin_x86_64/include
```

into

```
SDSL_BASE_DIR =  ~/matching_statistics/indexed_ms/sdsl-lite/build
INCLUDES      = $(SDSL_BASE_DIR)/include    
INCLUDES_PRIME =  /usr/lib  # or to the location of your standard libraries if in a different location
```
then build:

```
~$ cd ~/matching_statistics/indexed_ms/fast_ms
~$ make
~$ cd tests
~$ make  # this will run tests  
```

this should build the programs in `~/matching_statistics/indexed_ms/fast_ms/bin`

### Experiments

This section assumes that the installation procedure was completed successfully. Experiments are organized in folders under `~/matching_statistics/indexed_ms/tests/`:

 - `wl_tests/sae18`: tests related to rank optimizations for Weiner Link operations.
 - `wl_tests/sandbox_maxrep`: standbox tests related to rank optimizations for Weiner Link operations.
 - `parent_tests/sae18`: tests related to optimizations for parent operations.
 - `parent_tests/sandbox_selectatdist`: sandbox tests related to optimizations for parent operations.

Non-sandbox experiments, share input data hence they are generated only once under `~/matching_statistics/indexed_ms/tests/big_paper3/`. To generate the input data, type `make` under that directory.

To replicate an experiment, type `make` on its folder. This will generate `.csv` files with the results. We use `RStudio` (and the `tidyverse` meta-package) to generate plot and reports from the `.csv` files. Once, the experiments are completed, one can generate the reports under `~/matching_statistics/indexed_ms/tests/report`. 

# Contact
For suggestions or questions please write to Olgert Denas `gertidenas@gmail.com`
