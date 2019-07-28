# Data
Here we provide the inputs for the CellNOpt-ODE and the Dynamic-Feeder analysis.

These inputs consist of:
+ A **cnolist.RData** `R` object: Which contains the data obtained from perturbation experiments as described in [Hill et.al.](https://www.nature.com/articles/nmeth.3773)
+ A **model.RData** `R` object: Containing the Prior Knowledge Network (PKN) which is used to train to the data. The PKN was built and assembled as explained in [Razzaq et. al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006538)
+ A **database.RData** `R` object: Containing tha list of Kinase/Phosphatase-Substrate (K/P-S) knowledge and directed and signed Protein Interactions (PI's) from [OmniPath](http://omnipathdb.org/)
