setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
crit = 'AIC' # 'AIC'
crit = 'AIC'
library(doParallel) # start parallel computation

source('Functions/all_functions.R')
path.kdv <- ('pde_solver_data/kdv.mat')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

# kdv.mat <- readMat(path.kdv)
# kdv.t <- as.numeric(kdv.mat[['t']])
# kdv.x <- as.numeric(kdv.mat[['x']])
# kdv.usol <- as.matrix(Re(kdv.mat[['usol']]))
# dt = kdv.t[2]-kdv.t[1]
# dx = kdv.x[2]-kdv.x[1]
source("/nobackup/cfzh32/for_r/SG/kdv/kdv_solver.R")

candidate_library <- ASG_build_library(as.matrix(kdv.usol),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs

num <- round(10^(seq(2,4.8,0.2)))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
kdv_ada_lasso_pareto <- function(sam, ...){
  coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
  names(coeff_full) <- rhs
  coeff <- coeff_full[which(coeff_full!=0)]
  return(coeff)
}

system.time(kdv_ada_lasso_pareto_noiseless <- mclapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  for(j in 1:100){
    ada_lasso_pareto_out_list[[j]] <- kdv_ada_lasso_pareto(sample_index[[i]][[j]],crit=crit)
  }
  return(ada_lasso_pareto_out_list)
}, mc.cores = 6))
names(kdv_ada_lasso_pareto_noiseless) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
save(kdv_ada_lasso_pareto_noiseless, file=sprintf('Tests/Outputs/kdv_SG_increase_n_seed_%s_success_rate_pareto_%s.RData', seed,crit))

kdv_density_fun <- function(i, kdv_noiseless_data){
  kdv_noiseless <- kdv_noiseless_data[[i]]
  bernoulli_noiseless <- sapply(seq_along(kdv_noiseless), function(i){
    match_index <- match(names(kdv_noiseless[[i]]), kdv_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noiseless)
}
kdv_true_terms <- c("uu_{x}", "u_{xxx}")
kdv_density_noiseless_ada_lasso_pareto <- lapply(list(kdv_ada_lasso_pareto_noiseless), function(q) sapply(seq_along(q), kdv_density_fun, kdv_noiseless_data=q))


save(kdv_ada_lasso_pareto_noiseless, kdv_density_noiseless_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/kdv_SG_increase_n_seed_%s_success_rate_pareto_%s.RData', seed,crit))
