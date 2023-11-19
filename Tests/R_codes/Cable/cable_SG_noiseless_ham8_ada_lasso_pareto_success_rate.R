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

## Cable 
x = seq(-4,4,0.1)
t = seq(0,5,0.01)
init = exp(-(x)**2)
cable_out = pde_data(c(1,-1),c('u_{xx}','u'),x,t,init) # Cable 
cable_u <- t(as.matrix(cable_out[[1]]))
dx <- 0.1; dt <- 0.01
n <- length(t)
m <- length(x)

cable_true_terms <- c("u", "u_{xx}")

candidate_library <- ASG_build_library(as.matrix(cable_u),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs


# num <- 10^(seq(2,log10(20000),length=12))
num <- 10^(seq(2,4.4,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
# sample_index <- sapply(1:100, function(i) sample(nrow(R), n))
cable_ada_lasso_pareto <- function(sam, ...){
  coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
  names(coeff_full) <- rhs
  coeff <- coeff_full[which(coeff_full!=0)]
  return(coeff)
}

system.time(outs_ada_lasso_pareto_success_rate <- lapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  for(j in 1:100){
    out_temp <- tryCatch(out_temp <- cable_ada_lasso_pareto(sample_index[[i]][[j]]), 
             error = function(e) NULL)
    while(is.null(out_temp)){
      n_sam <- length(sample_index[[i]][[j]])
      sam_index <- sample(nrow(R), n_sam)
      out_temp <- tryCatch(out_temp <- cable_ada_lasso_pareto(sam_index), 
                           error = function(e) NULL)
    }
    ada_lasso_pareto_out_list[[j]] <- out_temp
  }
  return(ada_lasso_pareto_out_list)
}))
names(outs_ada_lasso_pareto_success_rate) <- sapply(1:length(num), function(j) paste0('n=',round(j)))

density_fun <- function(noise_coeff, true_terms){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(dens_noise)
}

cable_density_noiseless_ada_lasso <- sapply(seq_along(num), function(i) density_fun(outs_ada_lasso_pareto_success_rate[[i]], cable_true_terms))

save(outs_ada_lasso_pareto_success_rate,cable_density_noiseless_ada_lasso,
     file=sprintf('Tests/Outputs/cable_SG_noiseless_seed_%s_success_rate_ada_lasso_pareto.RData', seed))
