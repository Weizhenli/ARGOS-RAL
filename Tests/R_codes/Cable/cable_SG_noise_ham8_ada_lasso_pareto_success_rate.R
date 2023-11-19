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

density_noise_fun <- function(eta, u, num=100, seed=100, ...){
  # set.seed(seed)
  noise_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    un <- u + u_noise[[eta]][[i]]
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
system.time(cable_ada_lasso_pareto_noise <- parLapply(cl, seq_along(eta), density_noise_fun, u=cable_u, num=sample_size, seed=seed))
stopCluster(cl)
save(cable_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/cable_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))

cable_density_fun <- function(i, cable_noise_data){
  cable_noise <- cable_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(cable_noise), function(i){
    match_index <- match(names(cable_noise[[i]]), cable_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

# cable_bernoulli_all_noise_new <- lapply(cable_all_noise3, function(q) sapply(seq_along(q), cable_bernoulli_fun, cable_noise_data=q))
cable_density_noise_ada_lasso_pareto <- lapply(list(cable_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), cable_density_fun, cable_noise_data=q))

save(cable_ada_lasso_pareto_noise, cable_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/cable_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))
rm(list=ls())
