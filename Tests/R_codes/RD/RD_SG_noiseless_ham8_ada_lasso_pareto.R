setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate);library(doParallel)
seed = 10 # the same in python
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
num <- round(10^(seq(log10(200), log10(3000*50), length=12)))
num <- round(10^(seq(2.2, 5, 0.2)))
source('Functions/ada_lasso_pareto.R')
# sample data
num_xy = 5000; num_t=60
load(sprintf('pde_solver_data/RD_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
# identify 

fun_ada_lasso_pareto <- function(sam, ...){
  # df <- RD_noise_data[[12]] # noise=0
  # sam <- sample(1:nrow(df), j)
  # u_t
  RD_ut = ada_lasso_pareto_f(df[sam,-c(1,2)],df[sam,1], ...)
  # coef <- RD_ut$active_coeff
  # names(coef) <- names(df[,-c(1,2)])[which(w!=0)]
  coef_ut <- RD_ut$active_coeff
  # v_t
  RD_vt = ada_lasso_pareto_f(df[sam,-c(1,2)],df[sam,2], ...)
  # coef <- RD_vt$active_coeff
  # names(coef) <- names(df[,-c(1,2)])[which(w!=0)]
  coef_vt <- RD_vt$active_coeff
  return(list(ut=coef_ut, vt=coef_vt))
}

df <- RD_noise_data[[12]] # noise=0
s_time <- Sys.time()
# ada_lasso_pareto
system.time(
  RD_outs_ada_lasso_pareto_success_rate <- lapply(num, function(num_i){
    s_times = Sys.time()
    sam <- lapply(1:100, function(e) sample(nrow(df), num_i))
    out <- mclapply(1:100, mc.cores = 25, function(i) {
      out <-  fun_ada_lasso_pareto(sam[[i]], crit=crit); 
      #print(i); 
      return(out)
  })
    e_times = Sys.time()
    print(paste(num_i,difftime(e_times,s_times,unit='mins')))
    return(out)
  })
)
names(RD_outs_ada_lasso_pareto_success_rate) <- sapply(num, function(j) paste0('Number of points = ',j))
save(RD_outs_ada_lasso_pareto_success_rate, file=sprintf('Tests/Outputs/RD_SG_noiseless_seed_%s_ada_lasso_pareto.RData', seed)) # using hamilton

e_time <- Sys.time()
difftime(e_time,s_time)
