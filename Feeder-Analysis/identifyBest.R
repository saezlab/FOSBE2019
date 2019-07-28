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

mm = matrix(data = , nrow = 1, ncol = 3)
fits = c("1", "5", "10", "20")
pL = c("1", "2", "3", "4", "Inf")
penalty = c("2", "5", "10", "50", "100")

for(ii in 1:length(fits)){
  for(jj in 1:length(pL)){
    for(kk in 1:length(penalty)){
      mm = rbind(mm, t(as.matrix(c(fits[ii], pL[jj], penalty[kk]))))
    }
  }
}
mm = mm[-1, ]

idxBestAIC = NULL
aicBest = 100000
aicScores = c()
bicScores = c()
for(ii in 1:nrow(mm)){
  
  currFile = paste0("../Results/Cluster-Results/res_feeder_", mm[ii, 1], "_", mm[ii, 2], "_", mm[ii, 3], ".RData")
  if(file.exists(currFile)){
    
    load(file = currFile)
    
    score = aicCNO(model = res$`Integrated-Model`$model, cnolist = res$CNOList, opt_pars = res$Parameters)
    aicScores = c(aicScores, score$AIC)
    bicScores = c(bicScores, score$BIC)
    
    if(score$AIC < aicBest){
      
      aicBest = score$AIC
      idxBestAIC = ii
      
      save(res, file = "../Results/Best-Solutions/opt_pars_feeder.RData")
      
    }
    
  }
  
}

load(file = "../Results/Best-Solutions/opt_pars_feeder.RData")

source("../Public/map2cys.R")
attributes = map2cys(model = res$`Integrated-Model`$model, cnolist = res$CNOList, opt_pars = res$Parameters)
write.table(x = attributes$`Edge Attributes`, file = "../Results/Plots/Feeder-Model/edge_attributes_feeder.txt")
write.table(x = attributes$`Node Attributes`, file = "../Results/Plots/Feeder-Model/node_attributes_feeder.txt")
