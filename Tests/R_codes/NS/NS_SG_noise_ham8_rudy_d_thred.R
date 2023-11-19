setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate);
seed <- 10
set.seed(seed)
crit = 'AIC'

source('Functions/all_functions.R')


## load data and functions -------------------------------
num_xy = 5500; num_t = 65
load(sprintf('pde_solver_data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
num_xy = 5000; num_t = 60

sample_size <- 100

df_index <- lapply(seq_along(NS_noise_data), function(i){
  lapply(1:sample_size, function(x){
    sample(nrow(NS_noise_data[[i]]), num_xy*num_t)
  })
})



fun_rudy <- function(noise_level, df_index, lam=10^-5, d_tol=5, ...){
  # simulate data
  df <- NS_noise_data[[noise_level]][df_index, ]
  colnames(df)[2] <- '1'
  # rudy 
  w = TrainSTRidge(df[,-1],df[,1], lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df)[which(w!=0)+1]
  #beta2 <- coef.rudy
  return(coef.rudy)
}

## For noise density results --------------------
no_cores=12 # number of cpu cores that will bu used
# d_tols <- 1:10
d_tols <- c(0.2, 2, 10)
system.time(
  NS_noise_rudy_list <- mclapply(d_tols, function(tol){
    NS_noise_rudy <- mclapply(seq_along(NS_noise_data), function(k){
      # lapply(1:100, fun_rudy)
      rudy_list <- list()
      for(i in 1:sample_size){
        rudy_list[[i]] <- fun_rudy(k, df_index[[k]][[i]], d_tol=tol)
      }
      return(rudy_list)
    }, mc.cores = 6)
  }, mc.cores = 3)
)
save(NS_noise_rudy_list, file=sprintf('Tests/Outputs/NS_SG_noise_seed_%s_rudy_d_thred.RData', seed))

