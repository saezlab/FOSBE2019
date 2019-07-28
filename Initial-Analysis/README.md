# Initial-Analysis

Here we perform the training of the original PKN to data simply by executing the `initialAnalysis.R` script.

## `initialAnalysis.R`
+ [1-5]: Loading the necessary `R` packages.
+ [7-10]: Loading the PKN (`model.RData`) and the data (`cnolist.RData`) for the training.
+ [12-15]: Randomly setting the initial model parameter values (`k` and `tau` between `0` and `1` and fixed `n=3`).
+ [18-37]: Setting optimization pparamerter values.
+ [39]: Training the model to the data with the specified settings through an Enhanced Scatter Search (eSS) method.
+ [40]: Saving optimal model parameters identified by the optimization.
+ [42-46]: Assigning model parameters as network features for better visualization with Cytoscape.
+ [48-49]: Plotting and saving model simulations for the identified optimal set of parameters.
