[![DOI](https://zenodo.org/badge/383070693.svg)](https://zenodo.org/badge/latestdoi/383070693)

# hypocotyl
Python and C++ code to run a hypocotyl elongation model in A. thaliana.

Run main_figures.ipynb to get the simulation plots shown in the paper.

Run model_robustness.ipynb to get the supplementary figures related to 

Run hypocotyl_growth_model.ipynb to see how the model calculates and predicts the growth phenotypes of several genotypes.

File simulated_annealing_random.cpp holds the routine used to fit the model to the experimental dataset. 

Files functions_paper.* are header files for both C++ and Python codes. 

File heatmaps.cpp generates the data for Figure 6 in the manuscript. The simulated points used to generate the actual figure are in the results folder, but you can compile this file and check it for yourself.

The following files need to be compiled in order for the Notebooks to work:

1) File error.cpp yields the energy function associated with a given parameter set, so that we see the discrepancy between the fit and the data.

2) File sensitivity_analysis_global.cpp performs a sensitivity analysis of the model.

If you have any doubt or suggestion, feel free to contact me at: pcatalan [at] math [dot] uc3m [dor] com
