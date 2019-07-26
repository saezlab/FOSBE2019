library(CellNOptR)
library(MEIGOR)
library(CNORode)
library(doParallel)
library(readr)
library(infotheo)
library(igraph)

load(file = "../../Data/cnolist.RData")
load(file = "../../Data/model.RData")
load(file = "../../Data/database.RData")
load(file = "../../Analysis/simData.RData")

# model = preprocessing(data = cnolist, model = pknmodel, compression = TRUE, expansion = FALSE)

source("../../Public/computeMSE.R")
source("../../Public/computeMI.R")
source("../../Public/runDynamicFeeder.R")
source("../../Public/buildFeederObjectDynamic.R")
source("../../Public/identifyMisfitIndeces.R")
source("../../Public/map2cys.R")
source("../../Public/integrateLinks.R")
source("../../Public/preprocessingWeighted.R")

indeces <- identifyMisfitIndeces(cnolist = cnolist, model = model, simData = simData, mseThresh = 0.02)
object <- buildFeederObjectDynamic(model = model, cnolist = cnolist, database = database, indeces = indeces, pathLength = Inf)
integratedModel = integrateLinks(feederObject = object, cnolist = cnolist, compression = TRUE, expansion = FALSE, database = database)

plotModel(model = integratedModel$model, CNOlist = cnolist, indexIntegr = integratedModel$integLinksIdx)

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 15000;
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
paramsSSm$SSpenalty_fac=10
paramsSSm$SScontrolPenalty_fac=10

set.seed(5810)
# set initial parameters 
ode_parameters=createLBodeContPars(integratedModel$model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

res = runDynamicFeeder(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters, penFactor_k = 5, paramsSSm = paramsSSm)
save(res, file = "res_feeder.RData")

# plotLBodeFitness(cnolist = res$CNOList, model = res$Model, ode_parameters = res$Parameters, transfer_function = 4)

attributes = map2cys(model = res$Model, cnolist = res$CNOList, opt_pars = res$Parameters)

write.table(x = attributes$`Edge Attributes`, file = "edge_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = attributes$`Node Attributes`, file = "node_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
