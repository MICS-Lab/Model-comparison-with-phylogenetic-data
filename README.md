# Model choice for JAK2-mutated hematopoietic stem cells dynamics based on phylogenetic data

This repository contains all the scripts to build phylogenetic trees from genotype matrices and to run ABC-SMC and ABC-SMC model selection for the four patients presented. 

In the folder [data](data/), there are the genotype matrices for the four patients in the form of CSV files. 
The [input](input/) folder contains scripts to simulate the three models and to generate Lineage Through Time (LTT) plots and Clonal Fraction (CF) for the simulated tree. All these files are necessary to run the inference algorithm.

The scripts are made to run one step of the ABC-SMC algorithm for `N_particle` particles to approximate the sequential probability distributions. They are built to run in parallel across independent threads, improving computation time and flexibility.
To run one step of the model selection ABC-SMC for a single patient, refer to `(patient_ID)manual_model_selection.jl`. After each step run, it is possible to combine results across threads using `combine_smc.jl`, which also prints quantities of interest for use in the next step. 

To run a single model ABC-SMC, it is sufficient to change the model sampling in `smc_onestep.jl` as `m=1` for _VanEgeren et al._, `m=2` for _Williams et al._ and `m=3` for the adaptation to the construction of phylogenetic trees of _Hermange et al._. 

The script `run_allpatients.jl` is made to run model selection for the whole dataset for one step, then the results for one step across threads can be combined with `combine_allpatients.jl`.



