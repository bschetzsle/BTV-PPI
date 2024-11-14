#this script takes the posterior samples of partial correlation
#from the data, with separated stimulus effects,
#gets the median and makes the matrices symmetric,
#then saves it; it does not perform selection
setwd("~/BTV-PPI/Applied Results/separate_stimulus_effects")

make_symmetric <- function(X){
  P <- dim(X)[1]
  for(i in 1:(P-1) ){
    for(j in (i+1):P){
      X[i,j,] <- X[j,i,] <- (X[i,j,] + X[j,i,]) / 2
    }
  }
  return(X)
}

for(sub in 2:8){
  print(sprintf("Starting subject %s", sub))
  load(sprintf("subject_%s.RData", sub))
  N <- dim(result$par_cor_hat[[1]])[3]
  P <- dim(result$par_cor_hat[[1]])[1]
  K <- length(result$par_cor_hat) - 1
  n_samples <- dim(result$par_cor_hat[[1]])[4]
  
  for(k in 1:(K+1)){
    for(i in 1:n_samples){
      result$par_cor_hat[[k]][,,,i] <- result$par_cor_hat[[k]][,,,i] %>% 
        make_symmetric()
    }
  }

  par_cor_hat <- lapply(1:(K+1), function(k){
    apply(X = result$par_cor_hat[[k]], MARGIN = c(1,2,3), median)
  })
  
  save(par_cor_hat, file = sprintf("subject_%s_median.RData", sub))
}
