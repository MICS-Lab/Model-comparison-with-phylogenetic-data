# Model choice for JAK2-mutated hematopoietic stem cells dynamics based on phylogenetic data

This repository contains all the scripts to build phylogenetic trees from genotype matrices and to run ABC-SMC and ABC-SMC model selection for the four patients presented. 

In the folder [data](data/), there are the genotype matrices for the four patients in the form of CSV files. 
The [input](input/) folder contains scripts to simulate the three models and to generate Lineage Through Time (LTT) plots and Clonal Fraction (CF) for the simulated tree. All these files are necessary to run the inference algorithm.

To run the model selection ABC-SMC for a single patient, refer to `(patient_ID)manual_model_selection.jl`.



