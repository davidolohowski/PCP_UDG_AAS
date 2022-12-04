# Poisson Cluster Porcess Models for Detecting Ultra-Diffuse Galaxies 

This repository contains all source files, data, and figures required to reproduce the results contained in the paper "Poisson Cluster Porcess Models for Detecting Ultra-Diffuse Galaxies".

## Inference Algorithm

The MCMC inference algorithm for our models in the paper is contained in the file `src/BDMMH_supp.R` and `src/BDMMH_col.R`. The first one contains the algorithm for Model 1 while the second contains the one for Model 2.

## Model Fitting

For fitting the models to the data, run any one of the `pointing-ID.R` files. For example, run `v11acs.R` to fit our models to V11-ACS pointing. 

## Simulation

The simulations conducted in the paper are contained in the files `PCP_sim_XXX.R`. For example, run `PCP_sim_weak.R` to obtain simulations with weak UDG signal strengths.

## Data Analysis

For all post-model fitting data analysis, run any one of the `pointing-ID_da.R` file. For example, run `v11acs_da.R` to obtain all related posterior results and figures for the pointing V11-ACS.


