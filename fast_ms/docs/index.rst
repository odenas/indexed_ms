
The Snakemake Wrappers repository
=================================

Usage
-----

The general strategy is to include a wrapper into your workflow via the `wrapper <http://snakemake.readthedocs.io/en/latest/snakefiles/modularization.html?highlight=wrapper#wrappers>`_ directive, e.g.

.. code-block:: python

    rule all:
        input:
            s="inp.s",
            t="inp.t",
            fwd_cst="inp.s.fwd.stree",
            rev_cst="inp.s.rev.stree",
            maxrep="inp.s.rev.maxrep",
        output:
            "inp.t.ms"
        params:
            exec_path="<path to fast_ms>/bin/matching_stats_parallel.x",
            lca_parents=1,
            rank_fail=1,
            double_rank=1,
        threads:
            8
        wrapper:
            "file:" + str(Path(".").absolute().parent)
            "file:/" + "<path to fast_ms/wrappers/ms"


Here, Snakemake will automatically load the wrapper from the given path.

Each wrapper defines required software packages and versions. In combination
with the ``--use-conda`` flag of Snakemake, these will be deployed automatically.

.. toctree::
   :maxdepth: 3
   :glob:
   :hidden:
   :caption: Wrappers

   wrappers/*

.. toctree::
   :maxdepth: 3
   :glob:
   :hidden:
   :caption: Meta-Wrappers

   meta-wrappers/*
