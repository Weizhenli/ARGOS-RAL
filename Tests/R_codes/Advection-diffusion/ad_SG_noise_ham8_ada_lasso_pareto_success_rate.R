setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)

source_python('pde_solver_data/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

## advection-diffusion
x = seq(-10,10,0.1)
t = seq(0,10,0.01)
init = exp(-(x+2)**2)
ad_out = pde_data(c(1,-1),c('u_{xx}','u_{x}'),x,t,init) # advection-diffusion
ad_u <- t(as.matrix(ad_out[[1]]))
dx <- 0.1; dt <- 0.01
n <- length(t)
m <- length(x)

ad_true_terms <- c("u_{x}", "u_{xx}")

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(ad_u)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta, u, num=100, ...){
  # set.seed(seed)
  noise_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    # ad_un <- u + matrix(rnorm(m*n, 0, eta*sd(u)), nrow=m, ncol=n) # in linux it has problem
    ad_un <- u + u_noise[[eta]][[i]]
    candidate_library <- ASG_build_library(ad_un,dx,dt,noise=0)
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
# system.time(ad_all_noise <- parLapply(cl, 0:10*0.01, density_noise_fun, u=ad_u, num=sample_size, seed=seed))
system.time(ad_ada_lasso_pareto_noise <- parLapply(cl, seq_along(eta), density_noise_fun, u=ad_u, num=sample_size, seed=seed))
stopCluster(cl)
save(ad_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/ad_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))

ad_density_fun <- function(i, ad_noise_data){
  ad_noise <- ad_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(ad_noise), function(i){
    match_index <- match(names(ad_noise[[i]]), ad_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

# ad_bernoulli_all_noise_new <- lapply(ad_all_noise3, function(q) sapply(seq_along(q), ad_bernoulli_fun, ad_noise_data=q))
ad_density_noise_ada_lasso_pareto <- lapply(list(ad_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), ad_density_fun, ad_noise_data=q))

save(ad_ada_lasso_pareto_noise, ad_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/ad_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))
rm(list=ls())

