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
heat_ada_lasso_pareto <- function(sam, ...){
  coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
  names(coeff_full) <- rhs
  coeff <- coeff_full[which(coeff_full!=0)]
  return(coeff)
}

system.time(outs_ada_lasso_pareto_success_rate <- lapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  for(j in 1:100){
    temp_out <- tryCatch(temp_out <- heat_ada_lasso_pareto(sample_index[[i]][[j]]), 
                         error = function(e) 0)
    ada_lasso_pareto_out_list[[j]] <- temp_out
  }
  return(ada_lasso_pareto_out_list)
}))
names(outs_ada_lasso_pareto_success_rate) <- sapply(1:length(num), function(j) paste0('n=',round(j)))

density_fun <- function(noise_coeff, true_terms){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 1 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(dens_noise)
}

heat_density_noiseless_ada_lasso <- sapply(seq_along(num), function(i) density_fun(outs_ada_lasso_pareto_success_rate[[i]], heat_true_terms))

save(outs_ada_lasso_pareto_success_rate,heat_density_noiseless_ada_lasso,
     file=sprintf('Tests/Outputs/heat_SG_noiseless_seed_%s_success_rate_ada_lasso_pareto.RData', seed))
