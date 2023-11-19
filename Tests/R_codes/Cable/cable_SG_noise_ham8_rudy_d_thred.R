setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
crit = 'AIC'

source_python('pde_solver_data/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

cable_true_terms <- c("u", "u_{xx}")

# create noise 
snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

## Cable 
x = seq(-4,4,0.1)
t = seq(0,5,0.01)
init = exp(-(x)**2)
cable_out = pde_data(c(1,-1),c('u_{xx}','u'),x,t,init) # Cable 
cable_u <- t(as.matrix(cable_out[[1]]))
dx <- 0.1; dt <- 0.01
n <- length(t)
m <- length(x)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(cable_u)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta_i, u, num=100, lam=10^-5, d_tol=1, ...){
  # set.seed(seed)
  rudy_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    un <- u + u_noise[[eta_i]][[i]]
    candidate_library <- ASG_build_library(un,dx,dt,noise=0)
    Ut <- candidate_library$ut
    R <- candidate_library$Theta
    rhs <- candidate_library$rhs
    # rudy
    w = TrainSTRidge(R, Ut, lam, d_tol)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    rudy_list[[i]] <- coef.rudy
  }
  return(rudy_list)
}

sample_size <- 100
d_tols <- seq(0.2,2,0.2)
library(doParallel) # start parallel computation
no_cores <- 12 # number of cpu cores that will be used 
sample_size <- 100
# d_tols <- seq(0.2,2,0.2)
d_tols <- c(0.2, 2, 10)
library(doParallel) # start parallel computation
no_cores <- 22 # number of cpu cores that will be used 
system.time(
  cable_rudy_noise <- mclapply(d_tols, function(tol){
    mclapply(seq_along(eta), density_noise_fun, u=cable_u, num=sample_size, d_tol=tol, mc.cores=no_cores)
  }, mc.cores=3)
)
save(cable_rudy_noise, file=sprintf('Tests/Outputs/cable_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))

cable_density_fun <- function(cable_data){
  # cable_noise <- cable_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(cable_data), function(i){
    match_index <- match(names(cable_data[[i]]), cable_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

cable_density_noise_rudy <- lapply(cable_rudy_noise, function(k){
  sapply(k, cable_density_fun)
})

save(cable_rudy_noise, cable_density_noise_rudy,
     file=sprintf('Tests/Outputs/cable_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))
rm(list=ls())