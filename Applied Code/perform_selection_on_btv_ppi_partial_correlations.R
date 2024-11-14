library(tidyverse)

setwd("~/BTV-PPI/Applied Results")

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


for(subject in 1:8){
  print(sprintf("Performing selection on subject %s",subject))
  
  load(sprintf("subject_%s.RData", subject))
  
  pd <- apply(result$par_cor_hat, MARGIN = c(1,2,3), within_epsilon, epsilon = 0.05)
  
  decision_array <- fdr_select_pd(pd, fdr_target = 0.05)
  
  par_cor_btvppi <- apply(result$par_cor_hat, MARGIN = c(1,2,3), median)
  
  par_cor_btvppi_selection <- (par_cor_btvppi * decision_array) %>% 
    make_symmetric()
  
  result <- list(par_cor_btvppi = par_cor_btvppi,
                 par_cor_btvppi_selection = par_cor_btvppi_selection)
  
  save(result, file = sprintf("subject_%s_selection.RData", subject) )
}
