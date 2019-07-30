# Results

Here we provide all the results obtained from the analysis described in the other sections.

## Best-Solutions

+ `opt_pars_initial.RData` contains the model parameter after training the original PKN.
+ `simData_initial.RData` contains the model simulations for the optimal set of parameters identified after training the original PKN. This object is needed later as one of the inputs for the Dynamic-Feeder pipeline.
+ `opt_pars_feeder.RData` contains the model parameter after training the integrated model.

## Cluster-Results

Here we have provided the multiple integrated models that were evaluated by considering all the combinations of the following Dynamic-Feeder parameters.

+ Error thresholds at *5%*, *10%* and *20%* of the worse fitted measurements.
+ Maximal path length to search in the database for connecting the cues to the poorli fitted measurements. The path length values considered were *2*, *3*, *4* and *no path length limits*.
+ Penalty factor values over the newly integrated links at *5*, *10*, *50* and *100*.

These models were evaluated in a cluster infrastructure.

## Plots

### Initial-Model
Here we provide plots from the training of the initial PKN to data.
 
 + `Rplot - Initial Fit.pdf`: fit after training.
 + `Rplot - Initial Model.pdf`: initial PKN structure.
 + `initial_model_edge_attributes.txt`: inferred edge parameters (the *k* parameters)
 + `initial_model_node_attributes.txt`: inferred node parameters (the *tau* parameters)
 
### Feeder-Model
Here we provide plots from the training of the integrated network to data.
 
 + `Rplot - Feeder Fit.pdf`: fit after training.
 + `Rplot - Integrated Model.pdf`: integrated network structure.
 + `initial_model_edge_attributes.txt`: inferred edge parameters (the *k* parameters)
 + `initial_model_node_attributes.txt`: inferred node parameters (the *tau* parameters)
