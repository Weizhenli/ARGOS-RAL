rm(list=ls())
library(glmnet);library(reticulate)
library(doParallel)
seed=10
set.seed(seed)

# set path for NS-SG-sample.py
source_python('pde_solver_data/NS-SG-sample.py')
source('Functions/all_functions.R')

# simulate data with noise
fun_noise <- function(noise, num_xy, num_t, cores=0){
  df <- sam_ns(num_xy=num_xy, num_t=num_t, noise=noise, cores=cores)
  return(df)
}
num_xy = 5500; num_t = 65
snr_db_seq <- c(seq(0, 40, 2), Inf)
eta <- 10^(-snr_db_seq/20)

system.time(NS_noise_data <- mclapply(eta, fun_noise, cores=2, num_xy=num_xy, num_t=num_t, mc.cores=22)) # for noise plots
save(NS_noise_data,
     file=sprintf('pde_solver_data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
# user   system  elapsed 
# 91664.24  1363.84 17402.19 
# system.time(NS_noise_data <- lapply(0:5*0.01, fun_noise, cores=10, num_xy=5000, num_t=60)) # for conditional plots
# user     system    elapsed 
# 207833.253   1618.846  65315.815 
# names(NS_noise_data) <- sapply(0:5*0.01, function(x){paste0('noise=',x)})
names(NS_noise_data) <- sapply(snr_db_seq, function(x){paste0('SNR_dB=',x)})
#NS_names <- sam_ns(num_xy=100,num_t=1,noise=0.01,cores=0)
# save data
save(NS_noise_data,
     file=sprintf('pde_solver_data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))