setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
crit = 'AIC' # 'AIC'
crit = 'AIC'

source('Functions/all_functions.R')
path.kdv <- ('pde_solver_data/kdv.mat')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

kdv.mat <- readMat(path.kdv)
kdv.t <- as.numeric(kdv.mat[['t']])
kdv.x <- as.numeric(kdv.mat[['x']])
kdv.usol <- as.matrix(Re(kdv.mat[['usol']]))
dt = kdv.t[2]-kdv.t[1]
dx = kdv.x[2]-kdv.x[1]

n <- length(kdv.t) # 201
m <- length(kdv.x) # 512

kdv_true_terms <- c("uu_{x}", "u_{xxx}")

# create noise 
snr_db_seq <- c(seq(20, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(kdv.usol)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta, u, num=100, seed=100, ...){
  # set.seed(seed)
  noise_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    kdv_un <- u + matrix(rnorm(m*n, 0, eta*sd(u)), nrow=m, ncol=n) 
    candidate_library <- ASG_build_library(kdv_un,dx,dt,noise=0)
    Ut <- candidate_library$ut
    R <- candidate_library$Theta
    rhs <- candidate_library$rhs
    coeff_full <- ada_lasso_pareto_f(R, Ut, ...)[[1]]
    names(coeff_full) <- rhs
    coeff <- coeff_full[which(coeff_full!=0)]
    noise_list[[i]] <- coeff
  }
  return(noise_list)
}

sample_size <- 100
library(doParallel) # start parallel computation
no_cores <- 12 # number of cpu cores that will be used 
cl <- makeCluster(no_cores, type="FORK")
system.time(kdv_ada_lasso_pareto_noise <- parLapply(cl, eta, density_noise_fun, u=kdv.usol, num=sample_size, seed=seed, crit=crit))
stopCluster(cl)
save(kdv_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/kdv_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto_%s.RData', seed, sample_size,crit))

kdv_density_fun <- function(i, kdv_noise_data){
  kdv_noise <- kdv_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(kdv_noise), function(i){
    match_index <- match(names(kdv_noise[[i]]), kdv_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

kdv_density_noise_ada_lasso_pareto <- lapply(list(kdv_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), kdv_density_fun, kdv_noise_data=q))


save(kdv_ada_lasso_pareto_noise, kdv_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/kdv_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto_%s.RData', seed, sample_size,crit))
rm(list=ls())
