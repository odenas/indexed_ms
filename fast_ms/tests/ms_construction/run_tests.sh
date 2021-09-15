#!/bin/bash

# single threaded
snakemake -j1 --delete-all-output
snakemake -pr -j1

# multi threaded
snakemake -j1 --delete-all-output
snakemake -pr -j4
