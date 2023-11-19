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

## advection-diffusion
x = seq(-10,10,0.1)
t = seq(0,10,0.01)
init = exp(-(x+2)**2)
ad_out = pde_data(c(1,-1),c('u_{xx}','u_{x}'),x,t,init) # advection-diffusion
ad_u <- t(as.matrix(ad_out[[1]]))
dx <- 0.1; dt <- 0.01

ad_true_terms <- c("u_{x}", "u_{xx}")

candidate_library <- ASG_build_library(as.matrix(ad_u),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs

# num <- 10^(seq(2,log10(20000),length=12))
num <- 10^(seq(2,5.2,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
# sample_index <- sapply(1:100, function(i) sample(nrow(R), n))
ad_rudy <- function(sam, lam=10^-5, d_tol=1, ...){
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
        rudy_out_list[[j]] <- ad_rudy(sample_index[[i]][[j]], d_tol=tol)
      }
      return(rudy_out_list)
    }, mc.cores=8)
    names(outs) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
    return(outs)
  }, mc.cores=3)
)

ad_density_fun <- function(ad_data){
  # ad_noise <- ad_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(ad_data), function(i){
    match_index <- match(names(ad_data[[i]]), ad_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

success_rate_rudy_list <- lapply(outs_rudy_list, function(k){
  sapply(k, ad_density_fun)
})

save(outs_rudy_list, success_rate_rudy_list,
     file=sprintf('Tests/Outputs/more_test2/ad_SG_noiseless_seed_%s_success_rate_rudy_d_thred.RData', seed))
