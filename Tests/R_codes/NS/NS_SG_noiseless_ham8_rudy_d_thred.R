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

num_xy = 5500; num_t = 65
load(sprintf('Data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
num_xy = 5000; num_t = 60


fun_rudy <- function(j, sam, lam=10^-5, d_tol=5){
  # simulate data
  #df <- sam_ns(num_xy=5000,num_t=60,noise=noise)
  df <- NS_noise_data[[22]]
  # sam <- sample(1:nrow(df),j)
  colnames(df)[2] <- '1'
  # rudy 
  w = TrainSTRidge(df[sam,-1], df[sam,1], lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df)[which(w!=0)+1]
  #beta2 <- coef.rudy
  return(coef.rudy)
}
no_cores=15 # the number of cpu cores that will bu used 


num <- 10^(seq(log10(100), log10(60*5000), length=12))
num <- 10^(seq(2, 5.4, 0.2))
nrow_data <- nrow(NS_noise_data[[12]])
sam_all <- lapply(num, function(num_i){
  lapply(1:100, function(x){
    sample(1:nrow_data, num_i)
  })
})

# d_tols <- 1:10
d_tols <- c(0.2, 2, 10)
system.time(
  NS_outs_rudy_list <- lapply(d_tols, function(tol){
    NS_outs_rudy <- lapply(seq_along(num), function(num_i){
      mclapply(1:100, mc.cores = 25, function(i) fun_rudy(j=num[num_i], sam=sam_all[[num_i]][[i]] , d_tol=tol))
    })
    names(NS_outs_rudy) <- sapply(num, function(j) paste0('Number of points = ',j))
    return(NS_outs_rudy)
  })
)


save(NS_outs_rudy_list, file=sprintf('Tests/Outputs/NS_SG_noiseless_density_seed_%s_rudy_d_thred.RData', seed))

rm(list=ls())

