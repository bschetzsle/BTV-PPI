#This script fits the gPPI model to all 60 data simulations
library(tidyverse)

setwd("~/BTV-PPI/Simulation Results")

fit_one_row_gppi = function(data, region = 1, N, P, K){
  predictor_regions <- (1:P)[! 1:P %in% c(region)]
  filter_out <- sapply(1:K, function(k){sprintf("S%sY%s",k,region)})
  data <- select(data, -any_of(filter_out)) %>% 
    rename(response = sprintf("Y%s", region))
  
  fit <- lm(response ~ ., data = data)
  
  coefficients <- summary(fit)$coefficients[,1] * 
    ifelse(summary(fit)$coefficients[,4] < 0.05, 1, 0)
  sigma2 <- sum(fit$residuals^2)/(N - (P - 1)*K -1)
  
  result <- array(NA, dim = c(P, N))
  for(i in 1:P){
    if(i == region){
      result[i,] <- rep(1/sigma2, N)
    }else{
      temp_sum <- 
        lapply(1:K, function(k){
          data[sprintf("S%s",k)] * coefficients[sprintf("S%sY%s", k, i)]
         }) %>% Reduce("+", .)
      temp_sum <- temp_sum +
        rep(coefficients[sprintf("Y%s",i)], N)
      result[i,] <- t(temp_sum/sigma2)
    }
  }
  return(result)
}

fit_all_rows_gppi <- function(data, N, P, K){
  result <- 
    lapply(1:P, function(p){
      fit_one_row_gppi(data, region = p, N = N, P = P, K = K)
    })
  
  #Now it's time to aggregate the results into an Omega matrix
  Omega_hat <- array(NA, dim = c(P,P,N))
  for(i in 1:P){
    Omega_hat[i,,] <- result[[i]]
  }
  #Now it's time to derive the partial correlation estimates
  par_cor_hat <- array(NA, dim = c(P,P,N))
  for(i in 1:P){
    for(j in 1:P){
      par_cor_hat[i,j,] <- 
        -(Omega_hat[j,i,] + Omega_hat[i,j,])/
        sqrt(Omega_hat[i,i,] * Omega_hat[j,j,])/2
    }
  }
  return(list(Omega_hat = Omega_hat, par_cor_hat = par_cor_hat))
}


n_sims <- 60
N <- 1000
P <- 15
K <- 1

for(sim in 1:n_sims){
  print(sprintf("Starting simulation %s", sim))
  #start with the binary stimulus
  load(sprintf("./data/Simulation_%s_binary_stimulus.Rdata", sim))
  
  result <- fit_all_rows_gppi(data, N, P, K)
  
  save(result, file = sprintf("./gPPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  
  #now take that same structure but apply a dynamic stimulus
  load(sprintf("./data/Simulation_%s_dynamic_stimulus.Rdata", sim))
  result <- fit_all_rows_gppi(data, N, P, K)
  
  save(result, file = sprintf("./gPPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
}

