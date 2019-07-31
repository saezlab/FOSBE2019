# Initial-Analysis

Here we perform the training of the original PKN to data simply by executing the `initialAnalysis.R` script.

## `initialAnalysis.R`
+ [1-5]: Loading the necessary `R` packages.
+ [7-10]: Loading the PKN (`model.RData`) and the data (`cnolist.RData`) for the training.
+ [14-18]: Randomly setting the initial model parameter values (`k` and `tau` between `0` and `1` and fixed `n=3`).
+ [20-39]: Setting optimization parameter values.
+ [41]: Training the model to the data with the specified settings through an Enhanced Scatter Search (eSS) method and local DHC solver.
+ [42]: Saving optimal model parameters identified by the optimization.
+ [44-48]: Assigning model parameters as network features for better visualization with Cytoscape.
+ [50-51]: Plotting and saving model simulations for the identified optimal set of parameters.
