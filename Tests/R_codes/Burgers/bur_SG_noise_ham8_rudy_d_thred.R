setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
crit = 'AIC'

source('Functions/all_functions.R')
path.burgers <- ('Data/burgers.mat')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

burgers.mat <- readMat(path.burgers)
burgers.t <- as.numeric(burgers.mat[['t']])
burgers.x <- as.numeric(burgers.mat[['x']])
burgers.usol <- as.matrix(Re(burgers.mat[['u']]))
dt = burgers.t[2]-burgers.t[1]
dx = burgers.x[2]-burgers.x[1]

n <- length(burgers.t) # 101
m <- length(burgers.x) # 256

snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(burgers.usol)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta_i, u, num=100, lam=10^-5, d_tol=1, ...){
  # set.seed(seed)
  rudy_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    bur_un <- u + u_noise[[eta_i]][[i]]
    candidate_library <- ASG_build_library(bur_un,dx,dt,noise=0)
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
  bur_rudy_noise <- mclapply(d_tols, function(tol){
    mclapply(seq_along(eta), density_noise_fun, u=burgers.usol, num=sample_size, d_tol=tol, mc.cores=no_cores)
  }, mc.cores=3)
)
save(bur_rudy_noise, file=sprintf('Tests/Outputs/bur_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))

bur_true_terms <- c("uu_{x}", "u_{xx}")
bur_density_fun <- function(bur_data){
  # bur_noise <- bur_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(bur_data), function(i){
    match_index <- match(names(bur_data[[i]]), bur_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

bur_density_noise_rudy <- lapply(bur_rudy_noise, function(k){
  sapply(k, bur_density_fun)
})
  
save(bur_rudy_noise, bur_density_noise_rudy,
     file=sprintf('Tests/Outputs/bur_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))
rm(list=ls())