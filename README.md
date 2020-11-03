# Fast matching statistics in small space

Implementation and analysis of algorithms for the matching statistics problem. 

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
```

this should build the programs in `./fast_ms/bin` (files ending in `.x`). Run a program without arguments to get a description
of what it does and how to use.

To make it easy to run experiements and organize work we support [snakemake](https://snakemake.readthedocs.io/en/stable/index.html). To
get started, create a python virtual environment with Python 3. Either using conda

```
~$ conda create --prefix ./myenv python=3.7
~$ conda activate ./myenv
```
or virtualenv
```
~$ virtualenv myenv -p python3.7
~$ . myenv/bin/activate
```

With the virtual environment set up and active, we install the dependencies and make the documentation and tutorials

```
~$ pip install -r requirements.txt
~$ cd fast_ms/docs && make html
```

To view the documentation and the examples on workflows, point a browser to `fast_ms/docs/_build/html/index.html`.

With the environment active, you can also run tests:

```
cd tests
make
```


# Contact
For suggestions or questions please write to Olgert Denas `gertidenas@gmail.com`
