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
num_s <- prime_factorize(4e+5)
num_t_index <- sample(1:length(num_s),sample(1:3,1))
num_t <- num_s[num_t_index]
num_xy <- num_s[-num_t_index]
while(prod(num_t)>125 | prod(num_xy)>10000){
  num_t_index <- sample(1:length(num_s),sample(1:3,1))
  num_t <- num_s[num_t_index]
  num_xy <- num_s[-num_t_index]
}

## load data and functions -------------------------
num_xy = 5000; num_t=60
load(sprintf('pde_solver_data/RD_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
num_xy = 5000; num_t=30

data_w <- tryCatch(load('Tests/temp/RD_SG_noise_index.RData'), error = function(e) NULL)
if(is.null(data_w)){
  df_index <- lapply(seq_along(RD_noise_data), function(i){
    lapply(1:100, function(x){
      sample(nrow(RD_noise_data[[i]]), num_xy*num_t)
    })
  })
  save('Tests/temp/RD_SG_noise_index.RData')
}

fun_rudy <- function(noise, df_index, lam=10^-5, d_tol=1){
  #set.seed(100)
  #sam <- sample(nrow(df),j)
  # df <- RD_noise_data[[noise]]
  # sam <- 1:nrow(df)
  # sam <- sample(nrow(df), num)
  df <- RD_noise_data[[noise]][df_index, ]
  colnames(df)[2] <- '1'
  # u_t
  w = TrainSTRidge(df[,-c(1,2)],df[,1],lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df[,-c(1,2)])[which(w!=0)]
  beta2_ut <- coef.rudy
  # v_t
  w = TrainSTRidge(df[,-c(1,2)],df[,2],lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- names(df[,-c(1,2)])[which(w!=0)]
  beta2_vt <- coef.rudy
  return(list(ut=beta2_ut, vt=beta2_vt))
}

# no_cores=length(RD_noise_data) # number of cpu cores that will bu used
# d_tols <- seq(0.2,2,0.2)
d_tols <- c(0.2, 2, 10)
system.time(
  RD_noise_rudy_list <- lapply(d_tols, function(tol){
    RD_noise_rudy <- lapply(seq_along(RD_noise_data), function(k){
      mclapply(1:100, function(i) fun_rudy(k, df_index[[k]][[i]], d_tol=tol), mc.cores = 25)
      # rudy_list <- list()
      # for(i in 1:100){
      #   rudy_list[[i]] <- fun_rudy(k, df_index[[k]][[i]])
      # }
      # return(rudy_list)
    })
  })
)

save(RD_noise_rudy_list, file=sprintf('Tests/Outputs/RD_SG_noise_density_seed_%s_snr_rudy_d_thred.RData', seed))
rm(list=ls())