setwd("~/BTV-PPI/Applied Results")
library(tidyverse)

fir_convolve <- function(y, h){
  get_H <- function(h, dim1, dim2){
    dimh <- length(h)
    H <- array(0, dim=c(dim1, dim2))
    i <- 1
    while(i <= dim1 & i <= dim2){
      len <- min(dim1-i+1, dimh)
      H[i:(i+len-1),i] <- h[1:len]
      i <- i + 1
    }
    return(H)
  }
  
  n <- length(y)
  
  H <- get_H(h, n, n)
  
  return(H %*% y)
}

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
  
  result <- lapply(1:(K+1), function(i){array(0, dim = c(P, N, n_samples))})
  for(i in 1:P){
    if(i == region){
      result[[1]][i,,] <- t(1/fit$sigma2)
    }else{
      for(sample in 1:n_samples){
        result[[1]][i,,sample] <- -fit$beta[[sprintf("beta_Y%s",i)]][sample,2:(N+1)]/fit$sigma2[sample,]
        for(k in 1:K){
          result[[(k+1)]][i,,sample] <- -(data[[sprintf('S%s',k)]] *
            fit$beta[[sprintf("beta_S%sY%s",k,i)]][sample, 2:(N+1)])/fit$sigma2[sample,]
        }
      }
    }
  }
  return(result)
}

fit_all_rows <- function(data, N, P, K, n_samples = 500, n_burn = 50){
  result <- lapply(1:P, function(i){
    print(sprintf("Starting fit for region %s", i))
    fit_one_row(data, region = i, N = N, P = P, K = K, 
                n_samples = n_samples, n_burn = n_burn)
  })
  
  #Now it's time to aggregate the results into an Omega matrix (seperated by stimulus effects)
  Omega_hat <- lapply(1:(K+1), function(k){array(0, dim = c(P,P,N,n_samples))})
  for(k in 1:(K+1)){
    for(i in 1:P){
      Omega_hat[[k]][i,,,] <- result[[i]][[k]]
    }
  }

  #Now it's time to derive the partial correlation estimates
  rand_index <- lapply(1:P, function(i){sample(1:n_samples)})
  par_cor_hat <- lapply(1:(K+1), function(k){array(NA, dim = c(P,P,N,n_samples))})
  for(k in 1:(K+1)){
    for(sample in 1:n_samples){
      for(i in 1:P){
        for(j in 1:P){
          par_cor_hat[[k]][i,j,,sample] <- 
            -Omega_hat[[k]][i,j,,rand_index[[i]][sample]]/
            sqrt(Omega_hat[[1]][i,i,,rand_index[[i]][sample]] *
                   Omega_hat[[1]][j,j,,rand_index[[j]][sample]])
        }
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

hrf_function <- function(t){
  (t^(5)*exp(-t)/gamma(6) - t^(15)*exp(-t)/gamma(16)/6)
}

h <- seq(from = 0, by = 2, length.out = 30) %>% 
  hrf_function()
h <- h/sum(h)

###

K <- 4

for(subject in 1:8){
  load(sprintf("./processed_data/sep_learned_subject_%s.RData", subject))
  P <- dim(model_data$Y)[2]
  N <- dim(model_data$Y)[1]
  
  Y <- model_data$Y
  
  s <- vector(mode = "list", length = K)
  s[[1]] <- model_data$X_ind[,1]
  s[[2]] <- model_data$X_ind[,2]
  s[[3]] <- model_data$entropy_slow
  s[[4]] <- model_data$entropy_fast
  
  S <- lapply(1:K, function(k){fir_convolve(s[[k]], h)})
  
  
  data <- 
    data.frame(Y = Y, 
               S = Reduce("cbind", S), 
               SY = lapply(1:K, function(k){
                 lapply(1:P, function(p){
                   S[[k]] * Y[,p]
                 }) %>% Reduce("cbind",.)
               }) %>% Reduce("cbind",.)
    )
  
  colnames(data) <- c(sprintf("Y%s", 1:P), 
                      sapply(1:K, function(k){sprintf("S%s",k)}), 
                      lapply(1:K, function(k){
                        sapply(1:P, function(p){
                          sprintf("S%sY%s",k,p)
                        })
                      }) %>% unlist())
  
  n_samples <- 250
  n_burn <- 50
  result <- fit_all_rows(data, N = N, P = P, K = K, 
                         n_samples = n_samples, n_burn = n_burn)
  
  save(result, file = sprintf("~/BTV-PPI/Applied Results/separate_stimulus_effects/subject_%s.RData", subject))
}