library(tidyverse)
library(doParallel)

setwd("~/BTV-PPI/Simulation Results")

fit_one_row_bspline <- function(data, region = 1, N, P, K, df = 10){
  predictor_regions <- (1:P)[! 1:P %in% c(region)]
  filter_out <- sapply(1:K, function(k){sprintf("S%sY%s",k,region)}) %>% 
    c(sprintf("Y%s", region), .)
  X <- select(data, -any_of(filter_out))
  Y <- select(data, sprintf("Y%s", region)) %>% unlist()
  basis_matrix <- bs(1:N, intercept = TRUE, df = df, degree = 4)
  
  n_predictors <- dim(X)[2]
  
  #now you need to make all your predictors
  #you have 'df' basis functions; each predictor needs to be multiplied by each basis function, the kroenecker product
  B_X <- lapply(1:n_predictors, function(i){
    sapply(1:df, function(j){
      basis_matrix[,j] * X[,i]
    })
  }) %>%
    Reduce("cbind", .)
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x = B_X, y = Y, alpha = 1, intercept = FALSE,
                        lambda = exp(seq(2, -10, length.out=100)) )
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  
  fit <- glmnet(B_X, Y, alpha = 1, lambda = best_lambda, 
                intercept = FALSE)
  
  #Now I have to recover the elements of the row of Omega
  #C stores the weights for each basis function for each predictor
  #I have to start at coefficient 2 so that I don't include the intercept
  C <- coef(fit)[2:(df*n_predictors+1)] %>% 
    matrix(nrow = df, ncol = n_predictors, byrow = FALSE)
  #temp stores the varying coefficient for each regressor
  temp <- basis_matrix %*% C
  colnames(temp) <- colnames(X)
  
  result <- array(NA, dim = c(N,P))
  v <- (Y - predict(fit, newx = B_X))^2 %>% 
    sum()/(N-1)
  for(i in 1:P){
    if(i == region){
      result[,i] <- 1/v
    }else{
      temp_sum <- 
        lapply(1:K, function(k){
          data[sprintf('S%s',k)] *
            temp[,sprintf("S%sY%s",k,i)]
        }) %>% 
        Reduce("+",.)
      temp_sum <- -(temp_sum + temp[,sprintf("Y%s",i)])/v
      result[,i] <- unlist(temp_sum)
    }
  }
  
  return(result)
}

fit_all_rows_bspline <- function(data, N, P, K){
  cluster <- parallel::makeCluster(min(P, parallel::detectCores()-2) )
  doParallel::registerDoParallel(cl = cluster)
  
  result <- foreach(i = 1:P, .packages = c("dplyr","splines","glmnet"),
                    .export = c("fit_one_row_bspline")) %dopar% {
                      fit_one_row_bspline(data, region = i, N = N, P = P, 
                                          K = K, df = 10)
                    }
  stopCluster(cluster)
  
  #Now it's time to aggregate the results into an Omega matrix
  Omega_hat <- array(NA, dim = c(P,P,N))
  for(i in 1:P){
    Omega_hat[i,,] <- t(result[[i]])
  }
  #Now it's time to derive the partial correlation estimates
  par_cor_hat <- array(NA, dim = c(P,P,N))
  for(i in 1:P){
    for(j in 1:P){
      par_cor_hat[i,j,] <- -Omega_hat[i,j,] / 
        sqrt(Omega_hat[i,i,]*Omega_hat[j,j,])
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

N <- 1000
P <- 15
K <- 1

for(sim in 1:n_sims){
  print(sprintf("Starting simulation %s", sim))
  #start with the binary stimulus
  load(sprintf("./data/Simulation_%s_binary_stimulus.Rdata", sim))
  
  result <- fit_all_rows_bspline(data, N = N, P = P, K = K)
  
  save(result, file = sprintf("./splinePPI_results/Simulation_%s_Omega_binary_stimulus.Rdata", sim))
  
  #now take that same structure but apply a dynamic stimulus
  load(sprintf("./data/Simulation_%s_dynamic_stimulus.Rdata", sim))
  
  result <- fit_all_rows_bspline(data, N = N, P = P, K = K)
  
  save(result, file = sprintf("./splinePPI_results/Simulation_%s_Omega_dynamic_stimulus.Rdata", sim))
}
