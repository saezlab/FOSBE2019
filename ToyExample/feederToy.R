library(CellNOptR)
library(MEIGOR)
library(CNORode)
library(doParallel)
library(readr)
library(infotheo)
library(igraph)

# Sourcing all the Dynamic-Feeder functions needed
source("../Public/computeMSE.R")
source("../Public/computeMI.R")
source("../Public/runDynamicFeeder.R")
source("../Public/buildFeederObjectDynamic.R")
source("../Public/identifyMisfitIndices.R")
source("../Public/map2cys.R")
source("../Public/integrateLinks.R")
source("../Public/preprocessingWeighted.R")

# loading the toy example
data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")

model = ToyModel
cnolist = CNOlist(CNOlistToy)

# set initial parameters (here parameters 'k' and 'tau' are optimised and 'n' fixed to 3)
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 60;
paramsSSm$maxeval = Inf;
paramsSSm$atol=1e-6;
paramsSSm$reltol=1e-6;
paramsSSm$nan_fac=1000;
paramsSSm$dim_refset=30;
paramsSSm$n_diverse=1000;
paramsSSm$maxStepSize=Inf;
paramsSSm$maxNumSteps=10000;
paramsSSm$transfer_function = 4;

paramsSSm$lambda_tau=0.1
paramsSSm$lambda_k=0.01
paramsSSm$bootstrap=F
paramsSSm$SSpenalty_fac=0
paramsSSm$SScontrolPenalty_fac=0

## Training of the initial model
opt_pars=parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
simData = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars, transfer_function = 4)

## Loading interactions from Omnipath
database <- read_delim("omnipathDB.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
database <- as.matrix(database)

## Identifying the mis-fits (measurements with mse worse than 0.05) and interactions from the database which we want to integrate
indices = identifyMisfitIndices(cnolist = cnolist, model = model, simData = simData, mseThresh = 0.05)
feederObject = buildFeederObjectDynamic(model = model, cnolist = cnolist, indices = indices, database = database, pathLength = 2)
integratedModel = integrateLinks(feederObject = feederObject, cnolist = cnolist, compression = FALSE, expansion = FALSE, database = database)

## Plotting the integrated model by highlighting in purple the new added links to the PKN
plotModel(model = integratedModel$model, CNOlist = cnolist, indexIntegr = integratedModel$integLinksIdx)

## Optimizing the integrated model and plotting the improved fits
ode_parameters=createLBodeContPars(integratedModel$model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## penalty factor of 2
res = runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, paramsSSm = paramsSSm, penFactor_k = 2)
plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)

## penalty factor of 100000
res = runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, paramsSSm = paramsSSm, penFactor_k = 100000)
plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)
