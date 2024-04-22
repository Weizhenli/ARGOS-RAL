setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate);
seed <- 10
set.seed(seed)
crit = 'AIC'

source('Functions/all_functions.R')
## load data and functions -------------------------------
num_xy = 5500; num_t = 65
load(sprintf('Data/NS_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
num_xy = 5000; num_t = 60

sample_size <- 100

data_w <- tryCatch(load('Tests/temp/NS_SG_noise_index.RData'), error = function(e) NULL)
if(is.null(data_w)){
  df_index <- lapply(seq_along(NS_noise_data), function(i){
    lapply(1:100, function(x){
      sample(nrow(NS_noise_data[[i]]), num_xy*num_t)
    })
  })
  save(df_index,file='Tests/temp/NS_SG_noise_index.RData')
}


source('Functions/ada_lasso_pareto.R')

fun_our <- function(noise_level, df_index, ...){
  df <- NS_noise_data[[noise_level]][df_index, ]
  colnames(df)[2] <- '1'
  # our method
  model <- ada_lasso_pareto_f(df[,-1], df[,1], 1,...)
  beta <- model$active_coeff
  return(beta)
}

sample_size <- 100
# cl <- makeCluster(no_cores, type="FORK")
NS_noise_ada_lasso_pareto = list()
gamma_max <- 3
system.time(
  for(k in seq_along(NS_noise_data)){
    a = Sys.time()
    our_ada_lasso_pareto <- mclapply(1:sample_size, mc.cores=25, function(i) fun_our(k, df_index[[k]][[i]], crit=crit))
    NS_noise_ada_lasso_pareto[[k]] <- our_ada_lasso_pareto
    
    save(NS_noise_ada_lasso_pareto,
         file=sprintf('Tests/Outputs/NS_SG_noise_density_seed_%s_snr_ada_lasso_pareto.RData', seed))
    b = Sys.time()
    diff_time = difftime(b, a)
    print(paste(k, diff_time))
  }
)


