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
# sample data
num_xy = 5000; num_t=60
load(sprintf('pde_solver_data/RD_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
# identify 
df <- RD_noise_data[[length(RD_noise_data)]] # noise=0

nrow_data <- nrow(RD_noise_data[[12]])
sam_all <- lapply(num, function(num_i){
  lapply(1:100, function(x){
    sample(1:nrow_data, num_i)
  })
})

fun_rudy <- function(j, sam, lam=10^-5, d_tol=1){
  # sam <- sample(nrow(df), j)
  # u_t
  w = TrainSTRidge(df[sam,-c(1,2)], df[sam,1],lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df[,-c(1,2)])[which(w!=0)]
  beta2_ut <- coef.rudy
  # v_t
  w = TrainSTRidge(df[sam,-c(1,2)], df[sam,2],lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df[,-c(1,2)])[which(w!=0)]
  beta2_vt <- coef.rudy
  return(list(ut=beta2_ut, vt=beta2_vt))
}

s_time <- Sys.time()
# d_tols <- seq(0.2,2,0.2)
d_tols <- c(0.2, 2, 10)
# rudy
system.time(
  RD_outs_rudy_list <- lapply(d_tols, function(tol){
    RD_outs_rudy_success_rate <- lapply(seq_along(num), function(num_i){
      mclapply(1:100, function(i) fun_rudy(j=num[num_i], sam=sam_all[[num_i]][[i]], d_tol=tol), mc.cores = 25)
    })
    names(RD_outs_rudy_success_rate) <- sapply(num, function(j) paste0('Number of points = ',j))
    return(RD_outs_rudy_success_rate)
  })
)
save(RD_outs_rudy_list, file=sprintf('Tests/Outputs/RD_SG_noiseless_seed_%s_rudy_d_thred.RData', seed)) # using hamilton

e_time <- Sys.time()
difftime(e_time,s_time)
rm(list=ls())