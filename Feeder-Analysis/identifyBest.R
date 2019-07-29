library(CellNOptR)
library(MEIGOR)
library(CNORode)
library(doParallel)
library(readr)
library(infotheo)
library(igraph)

load(file = "../Results/Best-Solutions/opt_pars_initial.RData")
load(file = "../Data/cnolist.RData")
load(file = "../Data/model.RData")

source("../Public/aicCNO.R")
initial = aicCNO(model = model, cnolist = cnolist, opt_pars = opt_pars_initial)
initialAIC = initial$AIC
initialBIC = initial$BIC
rm(opt_pars_initial)

# mm = matrix(data = , nrow = 1, ncol = 3)
# # fits = c("1", "2", "5", "10", "20")
# fits = c("1", "5", "10", "20")
# pL = c("1", "2", "3", "4", "Inf")
# penalty = c("2", "5", "10", "50", "100")

error = c(0.05, 0.1, 0.2)
pL = c(2, 3, 4, Inf)
penalty = c(5, 10, 50, 100)

mm = matrix(data = , nrow = 1, ncol = 3)
for(ii in 1:length(error)){
  for(jj in 1:length(pL)){
    for(kk in 1:length(penalty)){
      toBind = t(as.matrix(c(error[ii], pL[jj], penalty[kk])))
      mm = rbind(mm, toBind)
    }
  }
}
mm = mm[-1, ]

idxBestAIC = NULL
aicBest = 100000
aicScores = c()
for(ii in 1:nrow(mm)){
  
  currFile = paste0("../Results/Cluster-Results/res_feeder_", mm[ii, 1], "_", mm[ii, 2], "_", mm[ii, 3], ".RData")
  if(file.exists(currFile)){
    
    load(file = currFile)
    
    score = aicCNO(model = res$`Integrated-Model`$model, cnolist = res$CNOList, opt_pars = res$Parameters)
    aicScores = c(aicScores, score$AIC)
    
    if(score$AIC < aicBest){
      
      aicBest = score$AIC
      idxBestAIC = ii
      
      save(res, file = "../Results/Best-Solutions/opt_pars_feeder.RData")
      
    }
    
  }
  
}

load(file = paste0("../Results/Cluster-Results/res_feeder_", mm[idxBestAIC, 1], "_", mm[idxBestAIC, 2], "_", mm[idxBestAIC, 3], ".RData"))

pdf("../Results/Plots/Feeder-Model/Rplot - Feeder Fit.pdf")
plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)
dev.off()

pdf("../Results/Plots/Feeder-Model/Rplot - Integrated Model.pdf")
plotModel(model = res$`Integrated-Model`$model, CNOlist = res$CNOList, indexIntegr = res$`Integrated-Model`$integLinksIdx)
dev.off()

save(res, file = "../Results/Best-Solutions/opt_pars_feeder.RData")
