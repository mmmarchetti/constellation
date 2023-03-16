# Exemplary Constellation Plots
This repository contains supplementary information for the paper "Finding constellations in chemical space through core analysis" in Frontiers in Chemistry by J. Jesús Naveja and José L. Medina-Franco.

### Introduction
The repository includes a Jupyter notebook that demonstrates the obtention of constellation plots for a library of Akt inhibitors and another of DNMT inhibitors.

The main goal of this work is to identify analog series or "constellations" in chemical space by identifying shared cores among compounds.

### Requirements
Prior to running the Jupyter notebook, you need to curate and process databases with the get_cores.py program, which is used to identify cores.

#### The notebook also requires the following libraries:

RDKit
Scikit-learn
Numpy
Pandas
Matplotlib

### Usage
You can run the Jupyter notebook to obtain constellation plots. The tsneFP function in the notebook is used to compute fingerprints for a library of molecules and apply t-SNE for dimensionality reduction. The function takes a pandas dataframe, a string matching the column name where the IDs of compounds is saved, a fingerprint type ("Morgan" or "MACCS"), perplexity and n_iter parameters for the t-SNE function, and a minimum size for the analog series to consider.

The output of the function is a new pandas dataframe containing information about cores, chemical space, and connections among cores.

Note that the final figures in the manuscript were produced in LaTeX using the data obtained by running the code in this notebook.
