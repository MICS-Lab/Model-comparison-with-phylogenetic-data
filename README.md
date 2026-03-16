# A Bayesian framework for comparing models of clonal dynamics from phylogenetic data, applied to JAK2V617F-driven Myeloproliferative Neoplasms

This repository contains all the scripts to build phylogenetic trees from genotype matrices and to run ABC-SMC and ABC-SMC model selection for the four patients presented. 

The data in the form of genotype matrices (.CSV files) should be put in the folder [data](data/). 
For data about patient P1 (ET2) and P2 (ET1) refer to:
Van Egeren, Debra, Javier Escabi, Maximilian Nguyen, et al. 2021. “Reconstructing the Lineage Histories and Differentiation Trajectories of Individual Cancer Cells in Myeloproliferative Neoplasms.” Cell Stem Cell 28 (3): 3. https://doi.org/10.1016/j.stem.2021.02.001.

For the data regarding patient P3 (PD5163) and P4 (PD7271): 
Williams, Nicholas, Joe Lee, Emily Mitchell, et al. 2022. “Life Histories of Myeloproliferative Neoplasms Inferred from Phylogenies.” Nature 602 (7895): 162–68. https://doi.org/10.1038/s41586-021-04312-6.

The [input](input/) folder contains scripts to simulate the three models and to generate Lineage Through Time (LTT) plots and Clonal Fraction (CF) for the simulated tree. All these files are necessary to run the inference algorithm.
In particular, the script `cellSimulation.jl` contains all the functions to build the exact and approximated genealogical tree from the model M3, and it has been used to perform all the numerical analysis to assess the quality of the approximation. 
The mathematical definition of the model and its original application can be found: 
Hermange, Gurvan, Alicia Rakotonirainy, Mahmoud Bentriou, et al. 2022. “Inferring the Initiation and Development of Myeloproliferative Neoplasms.” Proceedings of the National Academy of Sciences 119 (37): 37. https://doi.org/10.1073/pnas.2120374119.

The scripts are made to run one step of the ABC-SMC algorithm for `N_points` particles to approximate the sequential probability distributions. They are built to run in parallel across independent threads, improving computation time and flexibility.
To run one step of the model selection ABC-SMC for a single patient, refer to `(patient_ID)manual_model_selection.jl`. After each step run, it is possible to combine results across threads using `combine_smc.jl`, which also prints quantities of interest for use in the next step. 

To run a single model ABC-SMC, it is sufficient to change the model sampling in `smc_onestep.jl` as `m=1` for _VanEgeren et al._, `m=2` for _Williams et al._ and `m=3` for the adaptation to the construction of phylogenetic trees of _Hermange et al._. 

The script `run_allpatients.jl` is made to run model selection for the whole dataset for one step, then the results for one step across threads can be combined with `combine_allpatients.jl`.

The notebook `tree_simulation.ipynb` contains examples to simulate clonal dynamics and phylogenetic trees for the three models. It has been used to simulate and obtain synthetic datasets with specific parameters. 


