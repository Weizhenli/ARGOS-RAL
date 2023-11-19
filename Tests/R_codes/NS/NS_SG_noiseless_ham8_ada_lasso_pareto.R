setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate);
seed <- 10
set.seed(seed)
crit = 'AIC'

source('Functions/all_functions.R')


is_prime <- function(num){
  if(num<2){
    return(FALSE)
  }else if(num==2 | num==3){
    return(TRUE)
  }else{
    for(i in 2:sqrt(num)){
      if((num %% i)==0){
        return(FALSE)
      }
    }
    return(TRUE)
  }
}
#is_prime(125)
prime_factorize <- function(num){
  if(num<2){
    return(0)
  }
  if(is_prime(num)){
    return(0)
  }else{
    pri_fac_list <- numeric()
    while(TRUE){
      if(is_prime(num)){
        pri_fac_list <- c(pri_fac_list,num)
        break
      }
      for(i in 2:(num%/%2)){
        if(is_prime(i) & (num%%i)==0){
          pri_fac_list <- c(pri_fac_list,i)
          num <- num %/% i
          break
        }
      }
    }
  }
  return(pri_fac_list)
}
library(parallel)
all_var = ls()
num <- 10^(seq(log10(100), log10(60*5000), length=12))
num <- 10^(seq(2, 5.4, 0.2))
source('Functions/ada_lasso_pareto.R')
num_xy = 5500; num_t = 65
load(sprintf('pde_solver_data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
num_xy = 5000; num_t = 60
fun_ada_lasso_pareto <- function(j, ...){
  # simulate data
  df <- NS_noise_data[[22]]
  sam <- sample(1:nrow(df),j)
  colnames(df)[2] <- '1'
  # ada_lasso_pareto 
  w = ada_lasso_pareto_f(df[sam,-1],df[sam,1],1,...)
  coef <- w$active_coeff  # w[which(w!=0)]
  return(coef)
}
system.time(
  NS_outs_ada_lasso_pareto_success_rate <- lapply(num, function(num_i){
    mclapply(1:100, mc.cores = 25, function(i) fun_ada_lasso_pareto(j=num_i,crit=crit))
  })
)
names(NS_outs_ada_lasso_pareto_success_rate) <- sapply(num, function(j) paste0('Number of points = ',j))
save(NS_outs_ada_lasso_pareto_success_rate, file=sprintf('Tests/Outputs/NS_SG_noiseless_density_seed_%s_ada_lasso_pareto_%s.RData', seed,crit))



