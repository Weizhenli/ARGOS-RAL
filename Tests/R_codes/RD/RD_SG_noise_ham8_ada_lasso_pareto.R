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
  # save('Tests/temp/RD_SG_noise_index.RData')
}

source('Functions/ada_lasso_pareto.R')

fun_our <- function(noise, df_index, ...){
  df <- RD_noise_data[[noise]][df_index, ]
  colnames(df)[3] <- '1'
  # u_t
  model_ut <- ada_lasso_pareto_f(df[,-c(1,2)], df[,1], ...)
  coef_ut <- model_ut$active_coeff
  # v_t
  model_vt <- ada_lasso_pareto_f(df[,-c(1,2)], df[,2], ...)
  coef_vt <- model_vt$active_coeff
  return(list(ut = coef_ut, vt = coef_vt))
}

sample_size <- 100
RD_noise_ada_lasso_pareto = list()
system.time(
  for(k in 1:12){
    a = Sys.time()
    our_ada_lasso_pareto <- mclapply(1:sample_size, function(i) {
      out <- fun_our(k, df_index[[k]][[i]], crit=crit); 
      # print(i); 
      return(out)
    }, mc.cores=25)
    
    RD_noise_ada_lasso_pareto[[k]] <- our_ada_lasso_pareto

    save(RD_noise_ada_lasso_pareto,
         file=sprintf('Tests/Outputs/RD_SG_noise_density_seed_%s_snr_ada_lasso_pareto_%s.RData', seed,crit))
    b = Sys.time()
    diff_time = difftime(b, a, units='mins')
    print(paste(k, diff_time))
  }
)

