Not yet ready for public consumption.


# JS_SAA

In the spirit of reproducible research, this repository contains all the code necessary for running the experiments and generating graphs from the paper 
> "Gupta, Vishal and Kallus, Nathan, Data-Pooling in Stochastic Optimization (March 2019)."

Data for these experiments was taken from the [Kaggle Challenge](https://www.kaggle.com/c/rossmann-store-sales).

If you find this code or the paper useful, ***please consider citing it***.

## Overview
All of the source code for computing solutions by various methods can be found in genPurpose.jl, written in Julia.  The file babyNewsvendor.jl contains specialized functions for the specific case of a newsvendor with bernoulli Demand.   These two files are wrapped together in the single module JS defined in JS_SAA.main.jl

### Data Cleaning
The files:
  - clean_data.R
  - BinningData.ipynb
perform the requisite data-cleaning (in R) and discretization used for some experiments (in juyter notebook, Julia).   

### Experiments
The files:
  - testNewsvendorPlots.jl
  - syntheticDataHarness.jl
  - test_SyntheticData.jl

call workhorse functions from the module JS to simulate data and run experiments used for generating plots for exploration or for the paper.  Notice that since these algorithms are ``embarassingly parallel", substantive speed-ups can be achieved by multi-threading.  The file test_SyntheticData.jl does this in a naive way way by batching.  

The pairs of files: 
  - rollingHarness.jl
  - test_Rolling.jl
simiarly run a series of experiments on the historical Rossman data.  

The file 
  - backtestHarness.jl
is deprecated but sets up a test also to be run on the historical Rossman data.

### Plotting 
The folder **RPlots** contains functions used to generate plots for the paper using R.    

### Misc
**ArchiveExperiments** contain experiments and files not used in the final draft.  **Notebooks** contains jupyter notebooks mostly used for algorithmic development, debugging and exploration.    


## Licensing

This code is available under the MIT License.  
Copyright (c) 2019 Vishal Gupta
