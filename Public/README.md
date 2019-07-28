# Public
Here are provided all the functions of the Dynamic-Feeder pipeline.

+ `aicCNO.R` takes the prior knowledge network, the data and model parameters as an input in order to estimate the Akaike (AIC) score of that specific model.
+ `buildFeederObjectDynamic.R` estimates the possible mechanisms of interactions to be added to the PKN from a database of interactions for improving the fitting cost.
+ `computeMI.R` identifies poorly fitted measurements for specific experimental conditions based on Mutual Information values between simulations and data.
+ `computeMSE.R` identifies poorly fitted measurements for specific experimental conditions based on the Mean Squared Error values between simulations and data.
+ `identifyMisfitIndeces.R` identifies poorly fitted measurements for specific experimental conditions.
+ `integrateLinks.R` it adds the new possible links identified and integrates them to the original PKN.
+ `map2cys.R` maps model parameters as network features for better visualization with Cytoscape.
+ `preprocessingWeighted.R` preprocesses the original and integrated PKN based on weights assigned on each edge.
+ `runDynamicFeeder.R` it optimizes the integrated model and evaluates the effects of the integrated links.
