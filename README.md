# CAMEOS
## Wang Lab CUMC

This repository contains code associated with the CAMEOS project. Main optimization code associated with the paper is available in "main/".

Some example data is already populated in the folder but to run for arbitrary genes of interest you will need some of
the additional helper code available in prepare/.

Further documentation is available in doc/.

A script relating to our design of RBS sequences can be found in aux/.

For a readable python implementation/explanation of the first-order optimization step, see the iPython notebook in this directory.

Requirements to run example code:
- Unix environment
- HDF5, gzip tools (packages on ubuntu: hdf5-tools, zlib1g-dev ... probably already present for many users)
- hmmer v3+ (installed via hmmer ubuntu package or from source, should be in PATH)
- julia v1+
    - with modules: BioAlignments, Logging, StatsBase, JLD, Distributions, BioSymbols, NPZ, ArgParse (see doc/ for more info)
- python

Requirements to run code on own proteins:
- all of the above requirements
- large multiple sequence alignments
    - (possible ballpark: N / sqrt(L) > ~200 if N is number of sequences in MSA, L is length of protein of interest)
- to train MRFs: CCMPred (https://github.com/soedinglab/CCMpred) or Gremlin.
- tensorflow + GPU if you want to quickly run local energy calculations
    - for individual pairs, CPU should be fine, a bit slower
    - local energy calculations are not required, but may be informative about best regions to double-encode proteins.
- a few extra calculations to get all parameters associated with model (see walk-through in doc/).
