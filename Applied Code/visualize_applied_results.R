library(tidyverse)

hrf_function <- function(t){
  (t^(5)*exp(-t)/gamma(6) - t^(15)*exp(-t)/gamma(16)/6)
}

h <- seq(from = 0, by = 2, length.out = 30) %>% 
  hrf_function()
h <- h/sum(h)

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

plot_df = function(df)
{
  df %>% mutate(t = 1:dim(df)[1]) %>% reshape2::melt(id.vars = "t") %>% 
    ggplot() + geom_line(aes(x = t, y = value, color = variable))
}

get_quantiles <- function(samples){
  #this function gets the credible interval based on a 2D array where
  #the first dimension is the number of Samples and the second dimension is N
  #it returns a 3xN array, the first row being upper, second row being median
  #and third row being lower
  n_samples <- dim(samples)[1]
  N <- dim(samples)[2]
  
  sapply(1:N, function(t){
    samples[,t] %>% quantile(probs = c(0.025, 0.5, 0.975))
  })
}

plot_quantiles <- function(quantiles, truth = NA, color = "blue", title = ""){
  N <- dim(quantiles)[2]
  p1 <- 
    data.frame(lower = quantiles[1,], 
               median = quantiles[2,], 
               upper = quantiles[3,],
               t = 1:N,
               truth = truth) %>% 
    ggplot() +
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill=color, alpha = 0.3) +
    geom_line(aes(x = t, y = median), color = color, alpha = 0.3) +
    ggtitle(title)
  if(!is.na(truth)){
    p1 <- p1 + geom_line(aes(x = t, y = truth), color = "black")
  }
  p1
}

plot_median_pc_vs_entropy <- function(i, j, slow = TRUE, selection = TRUE){
  if(selection){
    data <- 
      data.frame(pc = result$par_cor_btvppi_selection[i,j,],
                 entropy_slow = S[[3]],
                 entropy_fast = S[[4]])[-c(1:15),]
  }else{
    data <- 
      data.frame(pc = result$par_cor_btvppi[i,j,],
                 entropy_slow = S[[3]],
                 entropy_fast = S[[4]])[-c(1:15),]
  }
  
  if(slow){
    lm_model <- lm(pc ~ entropy_slow, data = data)$coefficients
    ggplot() +
      geom_point(aes(x = entropy_slow, y = pc), data = data) +
      geom_abline(intercept = lm_model[1], slope = lm_model[2], color = "red")
  }else{
    lm_model <- lm(pc ~ entropy_fast, data = data)$coefficients
    ggplot() +
      geom_point(aes(x = entropy_fast, y = pc), data = data) +
      geom_abline(intercept = lm_model[1], slope = lm_model[2], color = "red")
  }
}



setwd("~/BTV-PPI/Applied Results")

subject <- 2
load(sprintf("~/BTV-PPI/Applied Results/processed_data/sep_learned_subject_%s.RData", subject))
K <- 4
N <- dim(model_data$Y)[1]

s <- vector(mode = "list", length = K)
s[[1]] <- model_data$X_ind[,1]
s[[2]] <- model_data$X_ind[,2]
s[[3]] <- model_data$entropy_slow
s[[4]] <- model_data$entropy_fast

S <- lapply(1:K, function(k){fir_convolve(s[[k]], h)})

load(sprintf("~/BTV-PPI/Applied Results/subject_%s_selection.RData", subject))

index = c(10,11)
plot_median_pc_vs_entropy(index[1], index[2], slow = FALSE, selection = TRUE) +
  coord_cartesian(ylim = c(-0.1,1)) +
  ggtitle(sprintf("Partial Correlation between %s and %s 
                  as a function of Entropy under a Fast Learning Rate",
                  model_data$roi_names[index[1]], model_data$roi_names[index[2]]))

plot_median_pc_vs_entropy(index[1], index[2], slow = TRUE, selection = TRUE) +
  coord_cartesian(ylim = c(-0.1,1)) +
  ggtitle(sprintf("Partial Correlation between %s and %s 
                  as a function of Entropy under a Slow Learning Rate",
                  model_data$roi_names[index[1]], model_data$roi_names[index[2]]))

index = c(10,11)
to_plot_df <- 
  lapply(1:8, function(subject){
    load(sprintf("~/BTV-PPI/Applied Results/subject_%s_selection.RData", subject))
    load(sprintf("~/BTV-PPI/Applied Results/processed_data/sep_learned_subject_%s.RData", subject))
    data.frame(subject = subject, 
               pc = result$par_cor_btvppi_selection[index[1], index[2],],
               entropy_slow = fir_convolve(model_data$entropy_slow, h),
               entropy_fast = fir_convolve(model_data$entropy_fast, h) )
  }) %>% Reduce("rbind", .)

ggplot(to_plot_df, aes(x = entropy_slow, y = pc)) +
  geom_point(aes(color = factor(subject) ), alpha = 0.2) +
  theme_minimal() +
  labs(color = "Subject") +
  geom_smooth(method = "gam", aes(color = factor(subject)), show.legend = FALSE) +
  geom_smooth(method = "gam", color = "black", show.legend = FALSE) +
  coord_cartesian(xlim = c(1.1,1.4)) +
  ggtitle(sprintf("Partial correlation between %s and %s 
                  as a function of Entropy under a slow learning rate", model_data$roi_names[index[1]], model_data$roi_names[index[2]]))



index = c(18,11)
to_plot_df <- 
  lapply(1:8, function(subject){
    load(sprintf("~/BTV-PPI/Applied Results/subject_%s_selection.RData", subject))
    load(sprintf("~/BTV-PPI/Applied Results/processed_data/sep_learned_subject_%s.RData", subject))
    data.frame(subject = subject, 
               pc = result$par_cor_btvppi_selection[index[1], index[2],],
               entropy_slow = fir_convolve(model_data$entropy_slow, h),
               entropy_fast = fir_convolve(model_data$entropy_fast, h))
  }) %>% Reduce("rbind", .)

ggplot(to_plot_df, aes(x = entropy_slow, y = pc)) +
  geom_point(aes(color = factor(subject) ), alpha = 0.2) +
  theme_minimal() +
  labs(color = "Subject") +
  geom_smooth(method = "gam", aes(color = factor(subject)), show.legend = FALSE) +
  geom_smooth(method = "gam", color = "black", show.legend = FALSE) +
  coord_cartesian(xlim = c(1.1,1.4)) +
  ggtitle(sprintf("Partial correlation between %s and %s 
                  as a function of Entropy under a slow learning rate", model_data$roi_names[index[1]], model_data$roi_names[index[2]]))



index = c(13,17)
to_plot_df <- 
  lapply(1:8, function(subject){
    load(sprintf("~/BTV-PPI/Applied Results/subject_%s_selection.RData", subject))
    load(sprintf("~/BTV-PPI/Applied Results/processed_data/sep_learned_subject_%s.RData", subject))
    data.frame(subject = subject, 
               pc = result$par_cor_btvppi_selection[index[1], index[2],],
               entropy_slow = fir_convolve(model_data$entropy_slow, h),
               entropy_fast = fir_convolve(model_data$entropy_fast, h))
  }) %>% Reduce("rbind", .)

ggplot(to_plot_df, aes(x = entropy_slow, y = pc)) +
  geom_point(aes(color = factor(subject) ), alpha = 0.2) +
  theme_minimal() +
  labs(color = "Subject") +
  geom_smooth(method = "gam", aes(color = factor(subject)), show.legend = FALSE) +
  geom_smooth(method = "gam", color = "black", show.legend = FALSE) +
  coord_cartesian(xlim = c(1.1,1.4)) +
  ggtitle(sprintf("Partial correlation between %s and %s 
                  as a function of Entropy under a slow learning rate", model_data$roi_names[index[1]], model_data$roi_names[index[2]]))


index = c(4,13)
to_plot_df <- 
  lapply(1:8, function(subject){
    load(sprintf("~/BTV-PPI/Applied Results/subject_%s_selection.RData", subject))
    load(sprintf("~/BTV-PPI/Applied Results/processed_data/sep_learned_subject_%s.RData", subject))
    data.frame(subject = subject, 
               pc = result$par_cor_btvppi_selection[index[1], index[2],],
               entropy_slow = fir_convolve(model_data$entropy_slow, h),
               entropy_fast = fir_convolve(model_data$entropy_fast, h))
  }) %>% Reduce("rbind", .)

ggplot(to_plot_df, aes(x = entropy_slow, y = pc)) +
  geom_point(aes(color = factor(subject) ), alpha = 0.2) +
  theme_minimal() +
  labs(color = "Subject") +
  geom_smooth(method = "gam", aes(color = factor(subject)), show.legend = FALSE) +
  geom_smooth(method = "gam", color = "black", show.legend = FALSE) +
  coord_cartesian(xlim = c(1.1,1.4)) +
  ggtitle(sprintf("Partial correlation between %s and %s 
                  as a function of Entropy under a slow learning rate", model_data$roi_names[index[1]], model_data$roi_names[index[2]]))




