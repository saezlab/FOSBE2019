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

# model = preprocessing(data = cnolist, model = pknmodel, compression = TRUE, expansion = FALSE)

source("aicCNO.R")
initial = aicCNO(model = model, cnolist = cnolist, opt_pars = opt_pars_initial)
initialAIC = initial$AIC
initialBIC = initial$BIC
rm(opt_pars_initial)

mm = matrix(data = , nrow = 1, ncol = 3)
# fits = c("1", "2", "5", "10", "20")
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

# save(mm, file = "Plots/indecesMatrix.RData")

idxBestAIC = NULL
idxBestBIC = NULL
idxBestFit = NULL
aicBest = 100000
bicBest = 100000
fitBest = 100000
aicScores = c()
bicScores = c()
for(ii in 1:nrow(mm)){
  
  # currFile = paste0("Results/misFit_", mm[ii, 1], "/pL_", mm[ii, 2], "/res_feeder_", mm[ii, 1], "_", mm[ii, 2], "_", mm[ii, 3], ".RData")
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
      
      # save(res, file = "Plots/res_aic.RData")
      
      # save(idxBestAIC, file = "Plots/idxBestAIC.RData")
      # 
      # pdf(file = paste0("Plots/aic_fit.pdf"), width = 22, height = 10)
      # plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)
      # dev.off()
      # 
      # pdf(file = paste0("Plots/aic_model.pdf"), width = 18, height = 14)
      # plotModel(model = res$`Integrated-Model`$model, CNOlist = res$CNOList, indexIntegr = res$`Integrated-Model`$integLinksIdx)
      # dev.off()
      
    }
    
    if(score$BIC < bicBest){
      
      bicBest = score$BIC
      idxBestBIC = ii
      
      # save(res, file = "Plots/res_bic.RData")
      
      # save(idxBestBIC, file = "Plots/idxBestBIC.RData")
      # 
      # pdf(file = paste0("Plots/bic_fit.pdf"), width = 22, height = 10)
      # plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)
      # dev.off()
      # 
      # pdf(file = paste0("Plots/bic_model.pdf"), width = 18, height = 14)
      # plotModel(model = res$`Integrated-Model`$model, CNOlist = res$CNOList, indexIntegr = res$`Integrated-Model`$integLinksIdx)
      # dev.off()
      
    }
    
    if(res$Parameters$ssm_results$fbest < fitBest){
      
      fitBest = res$Parameters$ssm_results$fbest
      idxBestFit = ii
      
      # save(res, file = "Plots/res_fit.RData")
      
      # save(idxBestFit, file = "Plots/idxBestFit.RData")
      # 
      # pdf(file = paste0("Plots/rss_fit.pdf"), width = 22, height = 10)
      # plotLBodeFitness(cnolist = res$CNOList, model = res$`Integrated-Model`$model, ode_parameters = res$Parameters, transfer_function = 4)
      # dev.off()
      # 
      # pdf(file = paste0("Plots/rss_model.pdf"), width = 18, height = 14)
      # plotModel(model = res$`Integrated-Model`$model, CNOlist = res$CNOList, indexIntegr = res$`Integrated-Model`$integLinksIdx)
      # dev.off()
      
    }
    
  }
  
}

