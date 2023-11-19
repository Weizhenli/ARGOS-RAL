setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
# crit = 'AIC' # 'AIC'
# crit = 'AIC'

source('Functions/all_functions.R')
path.kdv <- ('pde_solver_data/kdv.mat')
#pathname <- file.path(path,'kdv.mat')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

# kdv.mat <- readMat(path.kdv)
# kdv.t <- as.numeric(kdv.mat[['t']])
# kdv.x <- as.numeric(kdv.mat[['x']])
# kdv.usol <- as.matrix(Re(kdv.mat[['usol']]))
# dt = kdv.t[2]-kdv.t[1]
# dx = kdv.x[2]-kdv.x[1]

source("/nobackup/cfzh32/for_r/SG/kdv/kdv_solver.R")
n <- length(kdv.t) # 201
m <- length(kdv.x) # 512

# create noise 
snr_db_seq <- c(seq(20, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(kdv.usol)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta_i, u, d_tol, num=100, seed=100, ...){
  # set.seed(seed)
  rudy_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    kdv_un <- u + u_noise[[eta_i]][[i]]
    candidate_library <- ASG_build_library(kdv_un,dx,dt,noise=0)
    Ut <- candidate_library$ut
    R <- candidate_library$Theta
    rhs <- candidate_library$rhs
    # rudy
    w = TrainSTRidge(R, Ut, lam=10^-5, d_tol=d_tol)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    rudy_list[[i]] <- coef.rudy
  }
  return(rudy_list)
}

sample_size <- 100
# d_tols <- 1:10
# d_tols <- seq(10,100,10)
d_tols <- c(0.2, 2, 10)
library(doParallel) # start parallel computation
no_cores <- 12 # number of cpu cores that will be used 
system.time(
  kdv_rudy_noise <- mclapply(d_tols, function(tol){
    mclapply(seq_along(eta), density_noise_fun, u=kdv.usol, num=sample_size, d_tol=tol, mc.cores=no_cores)
  },mc.cores=3)
)
save(kdv_rudy_noise, file=sprintf('Tests/Outputs/kdv_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))

kdv_density_fun <- function(kdv_data){
  # kdv_noise <- kdv_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(kdv_data), function(i){
    match_index <- match(names(kdv_data[[i]]), kdv_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}
kdv_true_terms <- c("uu_{x}", "u_{xxx}")
kdv_density_noise_rudy_list <- lapply(kdv_rudy_noise, function(k){
  sapply(k, kdv_density_fun)
  # lapply(k, function(q) sapply(seq_along(q), kdv_density_fun, kdv_data=q))
})

save(kdv_rudy_noise, kdv_density_noise_rudy_list,
     file=sprintf('Tests/Outputs/kdv_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))
rm(list=ls())
