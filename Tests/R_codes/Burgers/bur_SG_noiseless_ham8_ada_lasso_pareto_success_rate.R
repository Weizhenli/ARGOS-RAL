setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
crit = 'AIC'
crit = 'AIC'
library(doParallel) # start parallel computation

source('Functions/all_functions.R')
path.burgers <- ('Data/burgers.mat')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

burgers.mat <- readMat(path.burgers)
burgers.t <- as.numeric(burgers.mat[['t']])
burgers.x <- as.numeric(burgers.mat[['x']])
burgers.usol <- as.matrix(Re(burgers.mat[['u']]))
dt = burgers.t[2]-burgers.t[1]
dx = burgers.x[2]-burgers.x[1]


candidate_library <- ASG_build_library(as.matrix(burgers.usol),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs

num <- 10^(seq(2,4.2,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
bur_ada_lasso_pareto <- function(sam, ...){
  coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
  names(coeff_full) <- rhs
  coeff <- coeff_full[which(coeff_full!=0)]
  return(coeff)
}

system.time(outs_ada_lasso_pareto_success_rate <- lapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  for(j in 1:100){
    ada_lasso_pareto_out_list[[j]] <- bur_ada_lasso_pareto(sample_index[[i]][[j]],crit=crit)
  }
  return(ada_lasso_pareto_out_list)
}))
names(outs_ada_lasso_pareto_success_rate) <- sapply(1:length(num), function(j) paste0('n=',round(j)))

debur_fun_bur <- function(noise_coeff, true_terms){
  debur_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(debur_noise)
}
true_term_bur = c("uu_{x}", "u_{xx}")
bur_density_noiseless_ada_lasso <- sapply(seq_along(num), function(i) debur_fun_bur(outs_ada_lasso_pareto_success_rate[[i]], true_term_bur))

save(outs_ada_lasso_pareto_success_rate,bur_density_noiseless_ada_lasso,
     file=sprintf('Tests/Outputs/bur_SG_noiseless_seed_%s_success_rate_ada_lasso_pareto_%s.RData', seed,crit))
