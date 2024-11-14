library(tidyverse)
library(doParallel)
library(sqldf)

setwd("~/BTV-PPI/Applied Results")

fir_convolve <- function(impulse_times, impulse_strengths, output_times){
  hrf_function <- function(t){
    if(t <= 0 || t >= 60){return(0)}
    (t^(5)*exp(-t)/gamma(6) - t^(15)*exp(-t)/gamma(16)/6)
  }
  
  #for each impulse, figure out what it's strength is at each of the output_times
  lapply(1:length(impulse_times), function(n){
    sapply(output_times, function(t2){
      impulse_strengths[n] * hrf_function(t2 - impulse_times[n])
    })
  }) %>% Reduce("+", .)
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
  
  #result <- array(NA, dim = c(P, N, n_samples))
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

###

K <- 2

for(subject in 1:8){
  load(sprintf("./processed_data/sep_learned_subject_%s.RData", subject))
  load(sprintf("./processed_data/behavior_subject_%s.Rdata", subject))
  
  #aggregate the data across the hemispheres
  Y <- sapply(1:9, function(p){
    (model_data$Y[,p] + model_data$Y[,(p+9)])/2
  })
  
  P <- dim(Y)[2]
  N <- dim(Y)[1]
  
  S <- vector(mode = "list", length = K)
  
  ##first predictor is entropy_difference
  temp1 <- data.frame(scan_times = model_data$scan_times)
  temp2 <- sqldf("select coalesce(lag(time) over (order by time), 0) as start_time,
        time as end_time, entropy_difference, self_self from behavior")
  
  entropy_difference <- sqldf("select scan_times, start_time, end_time, coalesce(entropy_difference, 0) as entropy_difference from temp1
        left join temp2 on (scan_times >= start_time and scan_times < end_time) ")$entropy_difference
  
  S[[1]] <- fir_convolve(impulse_times = model_data$scan_times,
                         impulse_strengths = entropy_difference,
                         output_times = model_data$scan_times)
  
  ##second predictor is self-self transitions
  S[[2]] <- fir_convolve(impulse_times = behavior$time,
                         impulse_strengths = behavior$self_self,
                         output_times = model_data$scan_times)
  
  rm(entropy_difference, temp1, temp2, behavior, model_data)
  
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
  
  save(result, file = sprintf("~/BTV-PPI/Applied Results/entropy_difference/subject_%s_entropy_difference_separate_stimulus_effects.RData", subject))
}