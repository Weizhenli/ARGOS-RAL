setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
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

# rudy
bur_rudy <- function(sam, lam=10^-5, d_tol=1, ...){
  w = TrainSTRidge(R[sam,], Ut[sam,], lam=lam, d_tol=d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- rhs[which(w!=0)]
  return(coef.rudy)
}

# d_tols <- seq(0.2,2,0.2)
d_tols <- c(0.2, 2, 10)
system.time(
  outs_rudy_list <- mclapply(d_tols, function(tol){
    outs <- lapply(1:length(num), function(i){
      rudy_out_list <- list()
      for(j in 1:100){
        rudy_out_list[[j]] <- bur_rudy(sample_index[[i]][[j]], d_tol=tol)
      }
      return(rudy_out_list)
    })
    names(outs) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
    return(outs)
  }, mc.cores=3)
)

bur_true_terms <- c("uu_{x}", "u_{xx}")
bur_density_fun <- function(bur_data){
  # bur_noise <- bur_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(bur_data), function(i){
    match_index <- match(names(bur_data[[i]]), bur_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

success_rate_rudy_list <- lapply(outs_rudy_list, function(k){
  sapply(k, bur_density_fun)
})

save(outs_rudy_list, success_rate_rudy_list,
     file=sprintf('Tests/Outputs/bur_SG_noiseless_seed_%s_success_rate_rudy_d_thred.RData', seed))
rm(list=ls())
