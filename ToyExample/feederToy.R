library(CellNOptR)
library(MEIGOR)
library(CNORode)
library(doParallel)
library(readr)
library(infotheo)
library(igraph)

# Sourcing all the functions needed
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
paramsSSm$maxtime = 900;
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

opt_pars=parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
simData = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars, transfer_function = 4)
