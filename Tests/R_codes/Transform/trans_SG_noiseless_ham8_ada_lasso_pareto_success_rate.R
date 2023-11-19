setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
library(doParallel) # start parallel computation

source_python('pde_solver_data/ode_int.py')
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

trans_true_terms <- c("u_{x}")

candidate_library <- ASG_build_library(as.matrix(trans_u),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs


# num <- 10^(seq(2,log10(20000),length=12))
num <- 10^(seq(2,5,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
# sample_index <- sapply(1:100, function(i) sample(nrow(R), n))
trans_ada_lasso_pareto <- function(sam, ...){
  coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
  names(coeff_full) <- rhs
  coeff <- coeff_full[which(coeff_full!=0)]
  return(coeff)
}

system.time(outs_ada_lasso_pareto_success_rate <- lapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  for(j in 1:100){
    ada_lasso_pareto_out_list[[j]] <- trans_ada_lasso_pareto(sample_index[[i]][[j]])
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

trans_density_noiseless_ada_lasso <- sapply(seq_along(num), function(i) density_fun(outs_ada_lasso_pareto_success_rate[[i]], trans_true_terms))

save(outs_ada_lasso_pareto_success_rate,trans_density_noiseless_ada_lasso,
     file=sprintf('Tests/Outputs/trans_SG_noiseless_seed_%s_success_rate_ada_lasso_pareto.RData', seed))
