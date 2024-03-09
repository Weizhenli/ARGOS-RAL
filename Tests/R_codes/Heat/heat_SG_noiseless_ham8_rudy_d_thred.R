setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
library(doParallel) # start parallel computation

source_python('PDE_solver/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

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

heat_true_terms <- c("u_{xx}")

candidate_library <- ASG_build_library(as.matrix(heat_u),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs

# num <- 10^(seq(2,log10(20000),length=12))
num <- 10^(seq(2,4.8,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
# sample_index <- sapply(1:100, function(i) sample(nrow(R), n))
heat_rudy <- function(sam, lam=10^-5, d_tol=1, ...){
  w = TrainSTRidge(R[sam,], Ut[sam,], lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- rhs[which(w!=0)]
  return(coef.rudy)
}

d_tols <- c(0.2, 2, 10)
system.time(
  outs_rudy_list <- mclapply(d_tols, function(tol){
    outs <- mclapply(1:length(num), function(i){
      rudy_out_list <- list()
      for(j in 1:100){
        rudy_out_list[[j]] <- heat_rudy(sample_index[[i]][[j]], d_tol=tol)
      }
      return(rudy_out_list)
    }, mc.cores=8)
    names(outs) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
    return(outs)
  }, mc.cores=3)
)

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

success_rate_rudy_list <- lapply(outs_rudy_list, function(k){
  sapply(k, heat_density_fun)
})

save(outs_rudy_list, success_rate_rudy_list,
     file=sprintf('Tests/Outputs/heat_SG_noiseless_seed_%s_success_rate_rudy_d_thred.RData', seed))
