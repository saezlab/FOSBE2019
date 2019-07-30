# FOSBE2019

Application of Dynamic-Feeder piepline over the HPN-DREAM dataset ([Hill et.al.]((https://www.nature.com/articles/nmeth.3773))).

## Organisation of the repository

+ [Data](https://github.com/saezlab/FOSBE2019/tree/master/Data): Here we provide the main inputs of our pipeline: a Prior Knowledge Network (PKN), a CNOlist and an Interaction Database object.
+ [Public](https://github.com/saezlab/FOSBE2019/tree/master/Public): Here we provide the main functions of the Dynamic-Feeder pipeline.
+ [Initial-Analysis](https://github.com/saezlab/FOSBE2019/tree/master/Initial-Analysis): Here we provide the script used for the training of the PKN to the data.
+ [Feeder-Analysis](https://github.com/saezlab/FOSBE2019/tree/master/Feeder-Analysis): Here we provide with the scripts used for the training of the integrated models and a script used for identifying the best model based on AIC score. By integrated model we are referring to the neworks which have new links identified from the Dynamic-Feeder pipeline and which are integrated to the PKN.
+ [Results](https://github.com/saezlab/FOSBE2019/tree/master/Results): Here we provide the object results together with the plots of the models we trained and the fits we obtained.

## Requirements

All network modelling steps were performed in [R](https://www.r-project.org/) v3.5.1 and visualised using [Cytoscape v3.4.](http://chianti.ucsd.edu/cytoscape-3.4.0/)

Other prerequisites include downloading and installing the following `R` package dependencies:

### Bioconductor

+ [CellNOptR](https://bioconductor.org/packages/release/bioc/html/CellNOptR.html)
+ [MEIGOR](https://www.bioconductor.org/packages/release/bioc/html/MEIGOR.html)
+ [CNORode](https://github.com/saezlab/CNORode)

### CRAN
+ [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
+ [readr](https://cran.r-project.org/web/packages/readr/index.html)
+ [infotheo](https://cran.r-project.org/web/packages/infotheo/index.html)
+ [igraph](https://cran.r-project.org/web/packages/igraph/index.html)

## How to run the pipeline

For better understanding we have compartmentalized the execution of the of the Dynamic-Feeder pipeline into two parts:

+ The initial training of the original PKN. We can train the original model simply by running `initialAnalysis.R` script [here](https://github.com/saezlab/FOSBE2019/blob/master/Initial-Analysis/initialAnalysis.R). Through this script we can identify the optimal set of parameters representing dynamic features of the elements in the PKN and we will generate the simulated behaviour of our system for these specific set of parameters. The latter will be used as an input for the Dynamic-Feeder pipeline.
+ We can then integrate new links to the original PKN and then do the training by running the `HPN_Feeder.R` script found [here](https://github.com/saezlab/FOSBE2019/blob/master/Feeder-Analysis/HPN_Feeder.R). This scrip will then identify the new set of optimal parameters for the integrated network.

## License

Distributed under the [GNU GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.html).
