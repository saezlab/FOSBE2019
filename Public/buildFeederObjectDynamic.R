#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E.Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

# This function estimates the possible mechanisms of interactions to be added to the PKN
# from a database of interactions for improving the fitting cost.

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             An indices object (i.e. as returned from computeMSE.R function) indicating poor fits
#             A database of interactions
# Optional:   A path length parameter for the maximal path length of additional interactions (pathLength = 4 by default)
#             A path length parameter for the maximal path length of additional interactions between the cues (cuesPathLength = 4 by default)
#             A loop length parameter for the inference of feedback mechanisms (loopLength = NULL by default)

buildFeederObjectDynamic <- function(model = model, cnolist = cnolist, indices = indices, database = database, pathLength = 3, cuesPathLength = NULL, loopLength = NULL){
  
  ##
  # Transform path lengs from based on interactions to based on nodes
  if((pathLength < 1) || !is.numeric(pathLength)){
    stop("Path length should be numeric and equal or greater to 1")
  } else {
    pathLength = pathLength + 1
  }
  if(!is.null(cuesPathLength)){
    cuesPathLength = cuesPathLength + 1
  }
  
  ##
  # Compacting cnolist
  if(class(cnolist)=="CNOlist"){
    cnolist = compatCNOlist(object = cnolist)
  }
  
  ##
  # Identifying the possible links
  modelSIF <- model2sif(model = model)
  if(length(indices$indices) > 0){
    indices <- indices$indices
  } else {
    indices <- NULL
  }
  
  df <- as.data.frame(x = database[, c(1, 3)])
  gg <- graph_from_data_frame(d = df, directed = TRUE)
  adj <- get.adjacency(gg)
  
  # all shortest paths connecting the measurements with the perturbed cues in the indices list
  sP_all <- list()
  if(!is.null(pathLength) && !is.null(indices)){
    
    for(ii in 1:length(indices)){
      
      measurement <- cnolist$namesSignals[indices[[ii]][1]]
      
      inhSet <- cnolist$namesCues[which(cnolist$valueCues[indices[[ii]][2], ]==1)]
      
      if(length(inhSet) > 0){
        
        for(jj in 1:length(inhSet)){
          
            if((inhSet[jj]%in%rownames(adj)) && (measurement%in%rownames(adj))){
              
              sP <- get.all.shortest.paths(graph = gg, from = which(rownames(adj)==inhSet[jj]), to = which(rownames(adj)==measurement))
              
              if(length(sP[[1]]) > 0){
                
                if(length(sP[[1]][[1]]) <= pathLength){
                  
                  for(kk in 1:length(sP[[1]])){
                    
                    sP_all[[length(sP_all)+1]] <- sP[[1]][[kk]]
                    
                  }
                  
                }
                
              }
              
            }
          
        }
        
      }
      
    }
    
  }
  
  # Now connecting the cues
  if(!is.null(cuesPathLength)){
    
    for(ii in 1:length(cnolist$namesCues)){
      
      for(jj in 1:length(cnolist$namesInhibitors)){
        
        if(ii != jj){
          
          if((cnolist$namesCues[ii]%in%rownames(adj)) && (cnolist$namesInhibitors[jj]%in%rownames(adj))){
            
            sP <- get.all.shortest.paths(graph = gg, from = which(rownames(adj)==cnolist$namesCues[ii]), to = which(rownames(adj)==cnolist$namesInhibitors[jj]))
            
            if(length(sP[[1]]) > 0){
              
              if(length(sP[[1]][[1]]) <= cuesPathLength){
                
                for(kk in 1:length(sP[[1]])){
                  
                  sP_all[[length(sP_all)+1]] <- sP[[1]][[kk]]
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  ##
  # now creating the feedInteractions list containing the signed interactions between the species
  feedInteractions <- list()
  for(ii in 1:length(sP_all)){
    
    feederMatrix <- matrix(data = , nrow = 1, ncol = 3)
    if(length(sP_all[[ii]]) > 1){
      
      for(jj in 1:(length(sP_all[[ii]])-1)){
        
        ss <- rownames(adj)[sP_all[[ii]][jj]]
        tt <- rownames(adj)[sP_all[[ii]][jj+1]]
        
        idx <- intersect(x = which(modelSIF[, 1]==ss), y = which(modelSIF[, 3]==tt))
        
        if(length(idx)==0){
          
          cc <- t(as.matrix(c(ss, database[intersect(which(database[, 1]==ss), which(database[, 3]==tt))[1], 2], tt)))
          
          feederMatrix <- unique(rbind(feederMatrix, cc))
          
        } else {
          
          cc <- t(as.matrix(c(ss, modelSIF[intersect(which(modelSIF[, 1]==ss), which(modelSIF[, 3]==tt))[1], 2], tt)))
          
          feederMatrix <- unique(rbind(feederMatrix, cc))
          
        }
        
      }
      
      if(nrow(feederMatrix)==2){
        
        feedInteractions[[length(feedInteractions)+1]] <- cc
        
      } else {
        
        feedInteractions[[length(feedInteractions)+1]] <- feederMatrix[-1, ]
        
      }
      
    }
    
  }
  
  ## 
  # identifying and adding loop interactions if specified
  if(!is.null(loopLength)){
    
    getLoops <- function(sif = sif, loopLength = loopLength){
      
      sif1 <- sif
      
      g <- graph_from_data_frame(d = as.data.frame(sif1[, c(1, 3)]), directed = TRUE)
      
      adj <- get.adjacency(g)
      
      adjacency <- get.adjacency(g)
      
      get_loops <- function(adj, paths, maxlen){
        maxlen <- maxlen - 1
        nxt_paths <- list()
        for(path in paths$paths){
          for(nxt in adj[[path[length(path)]]]){
            nxt_path <- c(path, nxt)
            if(path[1] == nxt & min(path) == nxt){
              paths$loops <- c(paths$loops, list(nxt_path))
            }else if(!(nxt %in% path)){
              nxt_paths <- c(nxt_paths, list(nxt_path))
            }
          }
        }
        paths$paths <- nxt_paths
        if(maxlen == 0){
          return(paths)
        }else{
          return(get_loops(adj, paths, maxlen))
        }
      }
      
      adj <- list()
      loops <- list()
      maxlen <- loopLength
      
      for(v in V(g)){
        adj[[as.numeric(v)]] <- neighbors(g, v)
      }
      
      for(start in seq(length(adj))){
        loops <- c(loops, get_loops(adj, list(paths = list(c(start)), 
                                              loops = list()), maxlen)$loops)
      }
      
      if(length(loops) > 0){
        
        ll <- list()
        for(i in 1:length(loops)){
          
          ll[[length(ll)+1]] <- rownames(adjacency)[loops[[i]][1:(length(loops[[i]])-1)]]
          
        }
        
        return(ll)
        
      }
      else{return(NULL)}
    }
    
    loopsAll <- getLoops(sif = database, loopLength = loopLength)
    
    idx <- c()
    for(ii in 1:length(loopsAll)){
      
      intSpecies <- intersect(x = loopsAll[[ii]], y = cnolist$namesSignals)
      
      if(length(intSpecies) > 0){
        
        idx <- c(idx, ii)
        
      }
      
    }
    
    if(length(idx)==0){
      
      print("There is no loop involving any of the measurements for this specific loop length..")
      
    } else {
      
      loopsAll <- loopsAll[idx]
      
      loopsAll <- unique(loopsAll)
      
      for(ii in 1:length(loopsAll)){
        
        if(loopsAll[[ii]][1]%in%model$namesSpecies){
          
          index <- which(loopsAll[[ii]]%in%cnolist$namesSignals)
          
          logicMatrix <- expand.grid(rep(list(0:1), length(index)))
          logicMatrix <- as.matrix(logicMatrix)
          
          for(jj in 1:nrow(logicMatrix)){
            
            feederMatrix <- matrix(data = , nrow = 1, ncol = 3)
            
            for(kk in 1:length(loopsAll[[ii]])){
              
              if(kk != length(loopsAll[[ii]])){
                
                ss <- loopsAll[[ii]][kk]
                tt <- loopsAll[[ii]][kk+1]
                
                if(kk%in%index){
                  
                  if(logicMatrix[jj, which(index==kk)]==0){ sign <- "-1"} else { sign <- "1"}
                  feederMatrix <- unique(rbind(feederMatrix, t(as.matrix(c(ss, sign, tt)))))
                  
                } else {
                  
                  idx1 <- which(database[, 1]==ss)
                  idx2 <- which(database[, 3]==tt)
                  
                  idx <- intersect(x = idx1, y = idx2)[1]
                  
                  sign <- database[idx, 2]
                  
                }
                
                feederMatrix <- unique(rbind(feederMatrix, t(as.matrix(c(ss, sign, tt)))))
                
              } else {
                
                ss <- loopsAll[[ii]][kk]
                tt <- loopsAll[[ii]][1]
                
                idx1 <- which(database[, 1]==ss)
                idx2 <- which(database[, 3]==tt)
                
                idx <- intersect(x = idx1, y = idx2)[1]
                
                sign <- database[idx, 2]
                
                feederMatrix <- unique(rbind(feederMatrix, t(as.matrix(c(ss, sign, tt)))))
                
              }
              
            }
            
            feedInteractions[[length(feedInteractions)+1]] <- feederMatrix[-1, ]
            
          }
          
        }
        
      }
      
    }
    
  }
  
  ## 
  # removing interactions which involve self-activating interactions
  if(length(feedInteractions)<=0){
    
    stop("No links to add for these settings..")
    
  } else {
    
    feedInteractions <- unique(feedInteractions)
    idx2rem <- c()
    for(ii in 1:length(feedInteractions)){
      if(class(feedInteractions[[ii]])=="matrix"){
        if(((feedInteractions[[ii]][1, 1]%in%model$namesSpecies)==FALSE) || (feedInteractions[[ii]][1, 1]==feedInteractions[[ii]][1, 3])){
          idx2rem <- c(idx2rem, ii)
        }
      } else {
        idx2rem <- c(idx2rem, ii)
      }
    }
    
    if(length(idx2rem)>0){
      feedInteractions <- feedInteractions[-idx2rem]
    }
    
    ##
    # removing those feedInteraction cases which are already present in the PKN
    idx2rem <- c()
    for(ii in 1:length(feedInteractions)){
      
      cnt = 0
      for(jj in 1:nrow(feedInteractions[[ii]])){
        
        idx1 = which(modelSIF[, 1]==feedInteractions[[ii]][jj, 1])
        idx2 = which(modelSIF[, 3]==feedInteractions[[ii]][jj, 3])
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx)>0){
          cnt <- cnt + 1
        }
        
      }
      
      if(cnt==nrow(feedInteractions[[ii]])){
        idx2rem <- c(idx2rem, ii)
        
      }
      
    }
    
    if(length(idx2rem)>0){
      feedInteractions <- feedInteractions[-idx2rem]
    }
    
    ##
    # removing mechanisms which might have incoming interactions in the stimuli's
    idx2rem = c()
    for(ii in 1:length(feedInteractions)){
      
      vv = unique(as.vector(feedInteractions[[ii]][, c(1, 3)]))[2:length(unique(as.vector(feedInteractions[[ii]][, c(1, 3)])))]
      
      int = intersect(x = vv, y = cnolist$namesStimuli)
      
      if(length(int) > 0){
        idx2rem = c(idx2rem, ii)
      }
      
    }
    
    if(length(idx2rem)>0){
      feedInteractions <- feedInteractions[-idx2rem]
    }
    
    ##
    # removing those feedInteractions cases which do not involve any cues
    idx2rem <- c()
    for(ii in 1:length(feedInteractions)){
      feedSpecies <- unique(c(feedInteractions[[ii]][, 1], feedInteractions[[ii]][, 3]))
      if(length(intersect(x = feedSpecies, y = cnolist$namesCues))==0){
        idx2rem <- c(idx2rem, ii)
      }
    }
    
    if(length(idx2rem)>0){
      feedInteractions <- feedInteractions[-idx2rem]
    }
    
    ##
    # keeping only those interactions which involve modt of the current present species in the PKN
    idx2keep <- c()
    dM <- matrix(data = , nrow = 1, ncol = 3)
    for(ii in 1:length(feedInteractions)){
      
      curr <- feedInteractions[[ii]]
      ss <- curr[1, 1]
      tt <- curr[nrow(curr), 3]
      currSpecies <- unique(c(curr[, 1], curr[, 3]))
      
      dM <- rbind(dM, t(as.matrix(c(ss, tt, length(intersect(x = currSpecies, y = model$namesSpecies)))))) 
      
    }
    
    dM <- dM[-1, ]
    
    if(is.null(nrow(dM))){
      dM = t(as.matrix(dM))
    }
    
    for(ii in 1:nrow(dM)){
      
      idx1 <- which(dM[, 1]==dM[ii, 1])
      idx2 <- which(dM[, 2]==dM[ii, 2])
      idx <- intersect(x = idx1, idx2)
      
      if(dM[ii, 3]==max(dM[idx, 3])){
        
        idx2keep <- c(idx2keep, ii)
        
      }
      
    }
    
    feedInteractions <- feedInteractions[idx2keep]
    
    if(length(feedInteractions)>0){
      
      ##
      # saving the object
      names(feedInteractions) <- paste0("Feed - ", 1:length(feedInteractions))
      
      object <- list()
      object[[length(object)+1]] <- modelSIF
      object[[length(object)+1]] <- feedInteractions
      
      names(object) <- c("Original PKN", "Feed mechanisms from database")
      
      return(object)
      
    } else {
      
      stop("No links to add for these settings..")
    }
    
  }
  
}