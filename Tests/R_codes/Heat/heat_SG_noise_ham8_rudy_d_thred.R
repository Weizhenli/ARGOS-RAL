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
  heat_rudy_noise <- mclapply(d_tols, function(tol){
    mclapply(seq_along(eta), density_noise_fun, u=heat_u, num=sample_size, d_tol=tol, mc.cores=no_cores)
  }, mc.cores=3)
)
save(heat_rudy_noise, file=sprintf('Tests/Outputs/heat_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))

heat_density_fun <- function(heat_data){
  # heat_noise <- heat_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(heat_data), function(i){
    match_index <- match(names(heat_data[[i]]), heat_true_terms)
    if(length(match_index) == 1 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

heat_density_noise_rudy <- lapply(heat_rudy_noise, function(k){
  sapply(k, heat_density_fun)
})

save(heat_rudy_noise, heat_density_noise_rudy,
     file=sprintf('Tests/Outputs/heat_SG_noise_seed_%s_samp_%s_snr_rudy_d_thred.RData', seed, sample_size))
rm(list=ls())