# FOSBE2019

Application of Dynamic-Feeder piepline over the [HPN-DREAM](https://www.nature.com/articles/nmeth.3773) dataset.

## Organisation of the repository

+ [Data](https://github.com/saezlab/FOSBE2019/tree/master/Data): Here we provide the main inputs of our pipeline: a Prioir Knowledge Network, a CNOlist and a Database object.
+ [Public](https://github.com/saezlab/FOSBE2019/tree/master/Public): Here we provide the main functions of the Dynamic-Feeder pipeline.
+ [Initial-Analysis](https://github.com/saezlab/FOSBE2019/tree/master/Initial-Analysis): Here we provide the script used for the inital training of the Prior Knowledge Network to the data.
+ [Feeder-Analysis](https://github.com/saezlab/FOSBE2019/tree/master/Feeder-Analysis): Here we provide with the scripts used for the training of the integrated models and a script used for identifying the best model based on AIC score.
+ [Results](https://github.com/saezlab/FOSBE2019/tree/master/Results): Here we provide the object results after the initial and Dynamic-Feeder training, together with plots of the models we train and the fits we obtain.

## Requirements

All network modelling steps were performed in [R](https://www.rstudio.com/) v3.5.1 and visualised using [Cytoscape](https://cytoscape.org/) v3.4.

Other prerequisites include downloading and installing the following `R` package dependencies:

+ [CellNOptR](https://bioconductor.org/packages/release/bioc/html/CellNOptR.html)
+ [MEIGOR](https://www.bioconductor.org/packages/release/bioc/html/MEIGOR.html)
+ [CNORode](https://github.com/saezlab/CNORode)
+ [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
+ [readr](https://cran.r-project.org/web/packages/readr/index.html)
+ [infotheo](https://cran.r-project.org/web/packages/infotheo/infotheo.pdf)
+ [igraph](https://igraph.org/r/)

The Analysis was run under MAC OS X El Captain machine.

## License

Distributed under the [GNU GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.html).

## How to run the pipeline

For better understanding we have compartmentalized the execution of the of the Dynamic-Feeder pipeline into two parts:

+ The initial training of the original PKN. We can train the original model simply by running `initialAnalysis.R` script [here](https://github.com/saezlab/FOSBE2019/blob/master/Initial-Analysis/initialAnalysis.R).
+ We can then integrate new links to the original PKN and then do the training by running the `HPN_Feeder.R` script found [here](https://github.com/saezlab/FOSBE2019/blob/master/Feeder-Analysis/HPN_Feeder.R)
