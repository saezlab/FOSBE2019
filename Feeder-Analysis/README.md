# Feeder-Analysis
Here we show the script which trains the integrated model to data and the one which identify the best model based on the AIC scores.

## `HPN_Feeder.R`
This is the script which trains the integrated model to the data. The integrated model was built by combining the original PKN with links.
Here the following Feeder parameters were used: path length limit `pL=Inf`; error threshold to evaluate `mseThresh=0.06291897` (corresponding to the value at 5% fo worst fit); edge penalty factor over the new integrated link equal to `penalty=5`.

+ [1-7]: Loading the necessary `R` packages.
+ [9-12]: Loading the PKN and the data; the database of interactions; and the simulations of the initial pkn with the best parameters that have been identified in the [Initial-Analysis](https://github.com/saezlab/FOSBE2019/tree/master/Initial-Analysis) directory.
+ [14-21]: Sourcing all the functions ncessary for running the Dynamic-Feeder pipeline.
+ [23]: Identifying the indeces of the measurements fitting worse than `mseThresh=0.06291897`. We can caluclate the MSE threshold value as: `sortedMSE = sort(x = as.vector(x = identifyMisfitIndeces(cnolist = cnolist, model = model, simData = simData)$mse), decreasing = TRUE); mseThresh = sortedMSE[round(0.05*length(sortedMSE))]`
+ [24]: Identifying all the paths connecting the cues to the measurements fitting worse than `mseThresh=06291897`.
+ [25]: Integrating those links to the original PKN.
+ [29-48]: Setting optimization parameters.
+ [51-55]: Setting initial model parameters randomly.
+ [57-59]: Optimizing model and saving results as an `R` object.
+ [61-64]: Mapping model parameters as network features and saving them as `.txt` file ready to by visualized in Cytoscape.
+ [66-72]: Saving Feeder fits and plot of the integrated model as pdf's.

## `identifyBest.R`
Script used to iterate accross all the combinations of Feeder parameters considered and then identifying and saving the one with the best AIC score.
