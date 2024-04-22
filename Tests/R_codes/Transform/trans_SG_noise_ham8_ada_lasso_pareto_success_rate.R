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

density_noise_fun <- function(eta_i, u, num=100, seed=100, ...){
  # set.seed(seed)
  noise_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    un <- u + u_noise[[eta_i]][[i]]
    candidate_library <- ASG_build_library(un,dx,dt,noise=0)
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
no_cores <- 22 # number of cpu cores that will be used 
cl <- makeCluster(no_cores, type="FORK")
system.time(trans_ada_lasso_pareto_noise <- parLapply(cl, seq_along(eta), density_noise_fun, u=trans_u, num=sample_size, seed=seed))
stopCluster(cl)
save(trans_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/trans_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))

trans_density_fun <- function(i, trans_noise_data){
  trans_noise <- trans_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(trans_noise), function(i){
    match_index <- match(names(trans_noise[[i]]), trans_true_terms)
    if(length(match_index) == 1 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

trans_density_noise_ada_lasso_pareto <- lapply(list(trans_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), trans_density_fun, trans_noise_data=q))

save(trans_ada_lasso_pareto_noise, trans_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/trans_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))
rm(list=ls())
