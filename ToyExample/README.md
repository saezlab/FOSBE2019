# ToyExample

Here we show how Dynamic-Feeder can be applied on a single execution script over a simple toy example ([ToyMMB](http://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html)).

## feederToy.R

Users can simply run the `feederToy.R` script to perform this simple case analysis.

+ [1-7]: Loading the required `R` packages.
+ [9-17]: Sourcing all the Dynamic-Feeder functions needed.
+ [19-24]: Loading the toy example.
+ [26-30]: Setting the initial parameters.
+ [32-51]: Setting optimization settings.
+ [53-55]: Training of the initial model.
+ [57-59]: Loading interactions from Omnipath.
+ [61-64]: Identifying the misfits and interactions from the database which we want to integrate.
+ [66-67]: Plotting the integrated model by highlighting in purple the newly added links to the PKN.
+ [69-81]: Optimizing the integrated model and plotting the improved fits.

Here we also show the effects of applying a low and a high penalty factor over the newly integrated links.
