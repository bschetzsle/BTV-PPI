library(tidyverse)

setwd("~/BTV-PPI/Simulation Results")

within_epsilon <- function(X, epsilon){
  return(mean(abs(X)<=epsilon))
}

fdr_select_pd <- function(pd, fdr_target){
  dims <- dim(pd)
  R <- dims[1]
  N <- dims[3]
  l_tri <- lower.tri(diag(R))
  d_array <- array(0, c(R,R,N))
  d_mat <- diag(R)
  thresh <- seq(0, 1, length.out = 10000)
  fdr_ach <- rep(0, N)
  for(tt in 1:N){
    p <- pd[,,tt][l_tri]
    fdr_star <- rep(0, length(thresh))
    
    for(i in 1:length(thresh)){
      num_select <- sum(p<thresh[i])
      if(num_select==0){
        fdr_star[i] <- 0
      } else{
        fdr_star[i] <- sum(p*(p<thresh[i]))/num_select
      }
      
    }
    ind <- which.min(abs(fdr_star-fdr_target))[1]
    final_thresh <- thresh[ind]
    fdr_ach[tt] <- fdr_star[ind]
    d <- 1*(p<final_thresh)
    d_mat[l_tri] <- d
    for(i in 2:R){
      for(j in 1:i){
        d_mat[j,i] <- d_mat[i,j]
      }
    }
    d_array[,,tt] <- d_mat
    d_mat <- diag(R)
    
  }
  
  return(d_array)
}

make_symmetric <- function(X){
  P <- dim(X)[1]
  for(i in 1:(P-1) ){
    for(j in (i+1):P){
      X[i,j,] <- X[j,i,] <- (X[i,j,] + X[j,i,]) / 2
    }
  }
  return(X)
}

load("Binary_stimulus_Omega.Rdata")
load("Dynamic_stimulus_Omega.Rdata")

#####
set.seed(123)
n_sims <- 60

for(sim in 1:n_sims){
  #You have btvPPI estimates of partial correlation
  #now you just need to do selection on them
  print(sprintf("Starting simulation %s", sim))
  
  #load the posterior samples for the binary stimulus
  load(sprintf("./btvPPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  
  #get partial correlation without selection
  par_cor_btvppi <- apply(result$par_cor_hat, MARGIN = c(1,2,3), median) %>% 
    make_symmetric()
  
  #get partial correlation with selection
  pd <- apply(result$par_cor_hat, MARGIN = c(1,2,3), within_epsilon, epsilon = 0.05)
  decision_array <- fdr_select_pd(pd, fdr_target = 0.05)
  par_cor_btvppi_selection <- (par_cor_btvppi * decision_array) %>% 
    make_symmetric()
  
  #now get the partial correlation for gPPI
  load(sprintf("./gPPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  par_cor_gppi <- result$par_cor_hat %>% make_symmetric()
  
  #now get the partial correlation for splinePPI
  load(sprintf("./splinePPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  par_cor_splineppi <- result$par_cor_hat %>% make_symmetric()
  
  result <- list(par_cor_btvppi = par_cor_btvppi,
                 par_cor_btvppi_selection = par_cor_btvppi_selection,
                 par_cor_gppi = par_cor_gppi,
                 par_cor_splineppi = par_cor_splineppi,
                 true_par_cor = par_cor)
  
  save(result, file = sprintf("./partial_correlations/Simulation_%s_par_cors_binary_stimulus.Rdata", sim) )
  
  #now repeat for the dynamic stimulus
  #load the posterior samples for the dynamic stimulus
  load(sprintf("./btvPPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
  
  #get partial correlation without selection
  par_cor_btvppi <- apply(result$par_cor_hat, MARGIN = c(1,2,3), median) %>% 
    make_symmetric()
  
  #get partial correlation with selection
  pd <- apply(result$par_cor_hat, MARGIN = c(1,2,3), within_epsilon, epsilon = 0.05)
  decision_array <- fdr_select_pd(pd, fdr_target = 0.05)
  par_cor_btvppi_selection <- (par_cor_btvppi * decision_array) %>% 
    make_symmetric()
  
  #now get the partial correlation for gPPI
  load(sprintf("./gPPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
  par_cor_gppi <- result$par_cor_hat %>% make_symmetric()
  
  #now get the partial correlation for splinePPI
  load(sprintf("./splinePPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
  par_cor_splineppi <- result$par_cor_hat %>% make_symmetric()
  
  result <- list(par_cor_btvppi = par_cor_btvppi,
                 par_cor_btvppi_selection = par_cor_btvppi_selection,
                 par_cor_gppi = par_cor_gppi,
                 par_cor_splineppi = par_cor_splineppi,
                 true_par_cor = par_cor_dynamic)
  
  save(result, file = sprintf("./partial_correlations/Simulation_%s_par_cors_dynamic_stimulus.Rdata", sim) )
}

