library(CellNOptR)
library(MEIGOR)
library(CNORode2017)
library(doParallel)
library(readr)

# loading the files
load(file = "../Data/cnolist.RData")
load(file = "../Data/database.RData")
load(file = "../Data/model.RData")

# set initial parameters (here parameters 'k' and 'tau' are optimised and 'n' fixed to 3)
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 3, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 7200;
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

opt_pars_initial=parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
save(opt_pars_initial, file = "../Results/Best-Solutions/opt_pars_initial.RData")

source(file = "../Public/map2cys.R")
attributes = map2cys(model = model, cnolist = cnolist, opt_pars = opt_pars_initial)

write.table(x = attributes$`Edge Attributes`, file = "../Results/Plots/Initial-Model/initial_model_edge_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = attributes$`Node Attributes`, file = "../Results/Plots/Initial-Model/initial_model_node_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

simData = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars_initial, transfer_function = 4)
save(simData, file = "../Results/Best-Solutions/simData_initial.RData")
