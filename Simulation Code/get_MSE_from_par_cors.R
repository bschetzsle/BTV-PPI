library(tidyverse)

make_symmetric <- function(X){
  P <- dim(X)[1]
  for(i in 1:(P-1) ){
    for(j in (i+1):P){
      X[i,j,] <- X[j,i,] <- (X[i,j,] + X[j,i,]) / 2
    }
  }
  return(X)
}

setwd("~/BTV-PPI/Simulation Results/partial_correlations")

situation_index <- list(c(1,13), c(4,7), c(1,2), c(13,14), c(10,11))
temp <- 
  lapply(1:60, function(simulation){
    load(sprintf("Simulation_%s_par_cors_binary_stimulus.Rdata",simulation))
    
    result$par_cor_gppi <- result$par_cor_gppi %>% make_symmetric()
    result$par_cor_btvppi <- result$par_cor_btvppi %>% make_symmetric()
    result$par_cor_splineppi <- result$par_cor_splineppi %>% make_symmetric()
    
    #calculate the SSE for current trial, for 5 situations, for 4 different models:
    #gPPI, splinePPI, btv-PPI, and btv-PPI with selection
    sapply(situation_index, function(index){
      rbind(
        (result$par_cor_gppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_splineppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_btvppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_btvppi_selection[index[1], index[2],] -
           result$true_par_cor[index[1], index[2],])^2 %>% sum()
      )
    }) %>% t()
  })

result <- array(NA,dim = c(5,4,60))
for(sim in 1:60){
  result[,,sim] <- temp[[sim]]
}

mean <- apply(result, MARGIN = c(1,2), mean)
sd <- apply(result, MARGIN = c(1,2), sd)

MSE_static <- list(mean = mean, sd = sd)
save(MSE_static, file = "MSE_static.Rdata")

temp <- 
  lapply(1:60, function(simulation){
    load(sprintf("Simulation_%s_par_cors_dynamic_stimulus.Rdata",simulation))
    
    result$par_cor_gppi <- result$par_cor_gppi %>% make_symmetric()
    result$par_cor_btvppi <- result$par_cor_btvppi %>% make_symmetric()
    result$par_cor_splineppi <- result$par_cor_splineppi %>% make_symmetric()
    
    #calculate the SSE for current trial, for 5 situations, for 4 different models:
    #gPPI, splinePPI, btv-PPI, and btv-PPI with selection
    sapply(situation_index, function(index){
      rbind(
        (result$par_cor_gppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_splineppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_btvppi[index[1], index[2],] - 
           result$true_par_cor[index[1], index[2],])^2 %>% sum(),
        (result$par_cor_btvppi_selection[index[1], index[2],] -
           result$true_par_cor[index[1], index[2],])^2 %>% sum()
      )
    }) %>% t()
  })

result <- array(NA,dim = c(5,4,60))
for(sim in 1:60){
  result[,,sim] <- temp[[sim]]
}

mean <- apply(result, MARGIN = c(1,2), mean)
sd <- apply(result, MARGIN = c(1,2), sd)

MSE_dynamic <- list(mean = mean, sd = sd)
save(MSE_dynamic, file = "MSE_dynamic.Rdata")
