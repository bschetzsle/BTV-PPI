library(tidyverse)

setwd("~/BTV-PPI/Simulation Results")

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

plot_quantiles <- function(quantiles, truth, color = "blue", title = ""){
  N <- dim(quantiles)[2]
  data.frame(lower = quantiles[1,], 
             median = quantiles[2,], 
             upper = quantiles[3,],
             t = 1:N,
             truth = truth) %>% 
    ggplot() +
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill=color, alpha = 0.3) +
    geom_line(aes(x = t, y = median), color = color, alpha = 0.3) +
    geom_line(aes(x = t, y = truth), color = "black") +
    ggtitle(title)
}

load("Sim_1_all_par_cors_static_stim.Rdata")

index = c(1,13)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Binary Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(4,7)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Binary Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(1,2)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Binary Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(13,14)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Binary Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(10,11)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Binary Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))


load("Sim_1_all_par_cors_dynamic_stim.Rdata")

index = c(1,13)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Dynamic Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(4,7)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Dynamic Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(1,2)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Dynamic Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(13,14)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Dynamic Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))

index = c(10,12)
data.frame(par_cor_btvppi_1 = result$par_cor_btvppi[index[1], index[2], ],
           par_cor_btvppi_2 = result$par_cor_btvppi[index[2], index[1], ],
           par_cor_btvppi_selection = result$par_cor_btvppi_selection[index[1], index[2], ],
           par_cor_gppi = -result$par_cor_gppi[index[1], index[2], ],
           true_par_cor = result$true_par_cor[index[1], index[2], ]) %>% 
  plot_df() +
  ggtitle(sprintf("Partial Correlation between regions %s and %s with a Dynamic Stimulus", index[1], index[2])) +
  coord_cartesian(ylim=c(-1,1))
