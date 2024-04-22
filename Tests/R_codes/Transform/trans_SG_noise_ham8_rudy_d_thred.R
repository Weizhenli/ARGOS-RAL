setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
crit = 'AIC'

source_python('PDE_solver/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

trans_true_terms <- c("u_{x}")

dx = 0.01; dt=0.01
x <- seq(-5,1,dx)
t <- seq(0,2,dt)
trans_u <- matrix(0,nrow=length(x),ncol=length(t))
for(i in 1:length(x)){
  for(j in 1:length(t)){
    trans_u[i,j] <- exp(-(x[i]+3*t[j])^2)
  }
}
n <- length(x); m <- length(t)
# create noise 
snr_db_seq <- c(seq(0, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(n*m, 0, noise*sd(trans_u)), nrow=n, ncol=m)
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
no_cores <- 11 # number of cpu cores that will be used 
system.time(
  trans_rudy_noise <- mclapply(d_tols, function(tol){
    mclapply(seq_along(eta), density_noise_fun, u=trans_u, num=sample_size, d_tol=tol, mc.cores=no_cores)
  }, mc.cores=3)
)
save(trans_rudy_noise, file=sprintf('Tests/Outputs/trans_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))


trans_density_fun <- function(trans_data){
  # trans_noise <- trans_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(trans_data), function(i){
    match_index <- match(names(trans_data[[i]]), trans_true_terms)
    if(length(match_index) == 1 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

trans_density_noise_rudy <- lapply(trans_rudy_noise, function(k){
  sapply(k, trans_density_fun)
})

save(trans_rudy_noise, trans_density_noise_rudy,
     file=sprintf('Tests/Outputs/trans_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))
rm(list=ls())