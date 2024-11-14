library(tidyverse)
library(doParallel)

setwd("~/BTV-PPI/Simulation Results")

fit_one_row <- function(data, region = 1, N, P, K, 
                        n_samples = 500, n_burn = 50){
  #this will return the samples of the a row of the precision matrix
  predictor_regions <- (1:P)[! 1:P %in% c(region)]
  filter_out <- sapply(1:K, function(k){sprintf("S%sY%s",k,region)})
  data <- select(data, -any_of(filter_out)) %>% 
    rename(response = sprintf("Y%s", region))
  
  fit <- shrinkTVP::shrinkTVP(response ~ ., 
                              data = data,
                              sv = TRUE,
                              niter = n_samples + n_burn, nburn = n_burn)
  
  result <- array(NA, dim = c(P, N, n_samples))
  for(i in 1:P){
    if(i == region){
      result[i,,] <- t(1/fit$sigma2)
    }else{
      for(sample in 1:n_samples){
        temp_sum <- 
          lapply(1:K, function(k){
            data[sprintf('S%s',k)] *
              fit$beta[[sprintf("beta_S%sY%s",k,i)]][sample, 2:(N+1)]
          }) %>% Reduce("+",.)
        temp_sum <- temp_sum + 
          fit$beta[[sprintf("beta_Y%s",i)]][sample,2:(N+1)]
        result[i,,sample] <- -temp_sum[,1]/fit$sigma2[sample,]
      }
    }
  }
  return(result)
}

fit_all_rows <- function(data, N, P, K, n_samples = 500, n_burn = 50){
  cluster <- parallel::makeCluster(min(P, parallel::detectCores()-2) )
  doParallel::registerDoParallel(cl = cluster)
  
  result <- foreach(i = 1:P, .packages = c("dplyr"),
                    .export = c("fit_one_row")) %dopar% {
                      fit_one_row(data, region = i, N = N, P = P, K = K, 
                                  n_samples = n_samples, n_burn = n_burn)
                    }
  stopCluster(cluster)
  
  #Now it's time to aggregate the results into an Omega matrix
  Omega_hat <- array(NA, dim = c(P,P,N,n_samples))
  for(i in 1:P){
    Omega_hat[i,,,] <- result[[i]]
  }
  #Now it's time to derive the partial correlation estimates
  rand_index <- lapply(1:P, function(i){sample(1:n_samples)})
  par_cor_hat <- array(NA, dim = c(P,P,N,n_samples))
  for(sample in 1:n_samples){
    for(i in 1:P){
      for(j in 1:P){
        par_cor_hat[i,j,,sample] <- 
          -Omega_hat[i,j,,rand_index[[i]][sample]]/
          sqrt(Omega_hat[i,i,,rand_index[[i]][sample]] *
                 Omega_hat[j,j,,rand_index[[j]][sample]])
      }
    }
  }
  return(list(Omega_hat = Omega_hat, par_cor_hat = par_cor_hat))
}

plot_df = function(df)
{
  df %>% mutate(t = 1:dim(df)[1]) %>% reshape2::melt(id.vars = "t") %>% 
    ggplot() + geom_line(aes(x = t, y = value, color = variable))
}


n_sims <- 60
set.seed(456)

n_samples <- 150
n_burn <- 50
n_thin <-  1

N <- 1000
P <- 15
K <- 1

for(sim in 1:n_sims){
  print(sprintf("Starting simulation %s", sim))
  #start with the binary stimulus
  load(sprintf("./data/Simulation_%s_binary_stimulus.Rdata", sim))
  
  result <- fit_all_rows(data, N = N, P = P, K = K, 
                         n_samples = n_samples, n_burn = n_burn)
  
  save(result, file = sprintf("./btvPPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  
  #now take that same structure but apply a dynamic stimulus
  load(sprintf("./data/Simulation_%s_dynamic_stimulus.Rdata", sim))
  
  result <- fit_all_rows(data, N = N, P = P, K = K, 
                         n_samples = n_samples, n_burn = n_burn)
  
  save(result, file = sprintf("./btvPPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
}
