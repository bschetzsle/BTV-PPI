library(tidyverse)

setwd("~/BTV-PPI/Simulation Results")

fir_convolve <- function(y, h){
  #this function convolves a stimulus vector y with the HRF h
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

plot_df = function(df){
  #this funtion plots the columns of a dataframe
  df %>% mutate(t = 1:dim(df)[1]) %>% reshape2::melt(id.vars = "t") %>% 
    ggplot() + geom_line(aes(x = t, y = value, color = variable))
}

test_positive_definite <- function(Omega){
  for(t in 1:N){
    if(min(eigen(Omega[,,t])$values) <= 0){
      return(FALSE)
    }
  }
  return(TRUE)
}

rcurve <- function(seed){
  #this function creates a random curve from basis splines
  set.seed(seed)
  basis_matrix <- splines::bs(1:N, knots = seq(1,N,length.out = 3), 
                              intercept = TRUE)
  X <- basis_matrix %*% matrix(runif(7, -1, 1), ncol = 1)
  X <- scales::rescale(X, to=c(max(min(X),-0.6),min(max(X),0.2) ) )
  return(X)
}

hrf_function <- function(t){
  (t^(5)*exp(-t)/gamma(6) - t^(15)*exp(-t)/gamma(16)/6)
}

h <- seq(from = 0, by = 2, length.out = 30) %>% 
  hrf_function()
h <- h/sum(h)

##First, generate the precision and covariance matrix from which the data will be simulated
set.seed(123)

N <- 1000
P <- 15
K <- 1

#determine s, both binary and dynamic
s <- vector(mode = "list", length = K)
s[[1]] <- rep(c(0,1), each = ceiling(N/4), length.out = N)
s_dynamic <- vector(mode = "list", length = K)
s_dynamic[[1]] <- 
  rnorm(N, mean = 0, sd = 1) %>% 
  cumsum() %>% 
  scales::rescale(to=c(0,1))

#convolve s with the hemodynamic response function
S <- vector(mode = "list", length = K)
S[[1]] <- fir_convolve(s[[1]], h)
S_dynamic <- vector(mode = "list", length = K)
S_dynamic[[1]] <- fir_convolve(s_dynamic[[1]], h)

#Construct the precision matrix directly, for both binary and dynamic S
Omega <- array(0, dim = c(P, P, N))
for(i in 1:P){
  Omega[i,i,] <- 1
}

#regions 1-3 are connected by a time-invariant physiological connection
for(i in 1:2){
  for(j in (i+1):3){
    Omega[i,j,] <- Omega[j,i,] <- runif(1,-0.3,0.3)
  }
}
test_positive_definite(Omega)

#regions 4-6 are connected by a time-invariant physiological connection
for(i in 4:5){
  for(j in (i+1):6){
    Omega[i,j,] <- Omega[j,i,] <- runif(1,-0.3,0.3)
  }
}
test_positive_definite(Omega)

#regions 7-9 are connected by a time-invariant physiological connection
for(i in 7:8){
  for(j in (i+1):9){
    Omega[i,j,] <- Omega[j,i,] <- runif(1,-0.3,0.3)
  }
}
test_positive_definite(Omega)

#regions 10-12 are connected by a varying physiological connection
for(i in 10:11){
  for(j in (i+1):12){
    Omega[i,j,] <- Omega[j,i,] <- rcurve(i*j)
  }
}
test_positive_definite(Omega)

#regions 13-15 are connected by a varying physiological connection
for(i in 13:14){
  for(j in (i+1):15){
    Omega[i,j,] <- Omega[j,i,] <- rcurve(i*j)
  }
}
test_positive_definite(Omega)

#regions 4 and 7 have only a PPI effect
Omega_dynamic <- Omega
Omega[4,7,] <- Omega[7,4,] <- S[[1]] * 0.3
Omega_dynamic[4,7,] <- Omega_dynamic[7,4,] <- S_dynamic[[1]] * 0.3
test_positive_definite(Omega)
test_positive_definite(Omega_dynamic)

#regions 10 and 11 have an additional PPI effect
Omega[10,11,] <- Omega[11,10,] <- rcurve(5) + S[[1]] * 0.3
Omega_dynamic[10,11,] <- Omega_dynamic[11,10,] <- rcurve(5) + S_dynamic[[1]] * 0.3
test_positive_definite(Omega)
test_positive_definite(Omega_dynamic)

#I think I have a positive-definite precision matrix!
par_cor <- par_cor_dynamic <- array(NA, dim = c(P, P, N))
for(t in 1:N){
  for(i in 1:P){
    for(j in 1:P){
      par_cor[i,j,t] <- -Omega[i,j,t] / sqrt(Omega[i,i,t] * Omega[j,j,t])
      par_cor_dynamic[i,j,t] <- 
        -Omega_dynamic[i,j,t] / sqrt(Omega_dynamic[i,i,t] * Omega_dynamic[j,j,t])
    }
  }
}

Sigma <- Sigma_dynamic <- array(NA, dim = c(P, P, N))
for(t in 1:N){
  Sigma[,,t] <- solve(Omega[,,t])
  Sigma_dynamic[,,t] <- solve(Omega_dynamic[,,t])
}

save(Omega = Omega, par_cor = par_cor, Sigma = Sigma, 
     file = "Binary_stimulus_Omega.Rdata")
save(Omega = Omega_dynamic, par_cor = par_cor_dynamic, Sigma = Sigma_dynamic,
     file = "Dynamic_stimulus_Omega.Rdata")

n_sims <- 60
for(sim in 1:n_sims){
  #first simulate data for the binary stimulus
  Y <- sapply(1:N, function(t){
    MASS::mvrnorm(1, mu = rep(0,P), Sigma = Sigma[,,t])
  })
  
  data <- cbind(t(Y),S[[1]],sapply(1:P, function(i){S[[1]]*Y[i,]}))
  colnames(data) <- c(sprintf("Y%s",1:P), "S1", sprintf("S1Y%s",1:P))
  data <- data.frame(data)
  
  save(data = data, 
       file = sprintf("./data/Simulation_%s_binary_stimulus.Rdata", sim))
  
  #now simulate data for the dynamic stimulus
  Y <- sapply(1:N, function(t){
    MASS::mvrnorm(1, mu = rep(0,P), Sigma = Sigma_dynamic[,,t])
  })
  
  data <- cbind(t(Y),S_dynamic[[1]],
                sapply(1:P, function(i){S_dynamic[[1]]*Y[i,]}))
  colnames(data) <- c(sprintf("Y%s",1:P), "S1", sprintf("S1Y%s",1:P))
  data <- data.frame(data)
  
  save(data = data, 
       file = sprintf("./data/Simulation_%s_dynamic_stimulus.Rdata", sim))
}
