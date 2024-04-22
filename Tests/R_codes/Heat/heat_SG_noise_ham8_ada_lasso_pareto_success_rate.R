setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)

source_python('PDE_solver/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

heat_true_terms <- c("u_{xx}")

# create noise 
snr_db_seq <- c(seq(40, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

# Heat
heat_sol <- function(x,t,k){
  L = max(x)
  u_sol <- matrix(0,nrow=length(x),ncol=length(t))
  for(i in 1:length(x)){
    for(j in 1:length(t)){
      u_sol[i,j] <- 6*sin(pi*x[i]/L)*exp(-k*(pi/L)^2*t[j])
    }
  }
  return(u_sol)
}
dx=0.01;dt=0.01
x <- seq(dx,5+dx,dx)
t <- seq(0,1.5,dt)
k <- 10
m <- length(x); n <- length(t)

# heat_u <- heat_sol(x,t,k)
init <- as.vector(heat_sol(x,0,2))
heat_out <- pde_data(list(1),list('u_{xx}'),x,t,init) 
heat_u <- t(as.matrix(heat_out[[1]]))

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(heat_u)), nrow=m, ncol=n)
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
system.time(heat_ada_lasso_pareto_noise <- parLapply(cl, seq_along(eta), density_noise_fun, u=heat_u, num=sample_size, seed=seed))
stopCluster(cl)
save(heat_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/heat_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))

heat_density_fun <- function(i, heat_noise_data){
  heat_noise <- heat_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(heat_noise), function(i){
    match_index <- match(names(heat_noise[[i]]), heat_true_terms)
    if(length(match_index) == 1 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

# heat_bernoulli_all_noise_new <- lapply(heat_all_noise3, function(q) sapply(seq_along(q), heat_bernoulli_fun, heat_noise_data=q))
heat_density_noise_ada_lasso_pareto <- lapply(list(heat_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), heat_density_fun, heat_noise_data=q))

save(heat_ada_lasso_pareto_noise, heat_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/heat_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto.RData', seed, sample_size))
rm(list=ls())
