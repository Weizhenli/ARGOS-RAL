setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal);library(caret)
seed = 10
set.seed(seed)
# crit = 'AIC' # 'AIC'
library(doParallel) # start parallel computation

source('Functions/all_functions.R')
path.kdv <- ('pde_solver_data/kdv.mat')
#pathname <- file.path(path,'kdv.mat')
cv <- import('cv2')
source('/nobackup/cfzh32/for_r/functions/MSA_lasso_pareto.R')

# kdv.mat <- readMat(path.kdv)
# kdv.t <- as.numeric(kdv.mat[['t']])
# kdv.x <- as.numeric(kdv.mat[['x']])
# kdv.usol <- as.matrix(Re(kdv.mat[['usol']]))
# dt = kdv.t[2]-kdv.t[1]
# dx = kdv.x[2]-kdv.x[1]

source("/nobackup/cfzh32/for_r/SG/kdv/kdv_solver.R")

candidate_library <- ASG_build_library(as.matrix(kdv.usol),dx,dt,noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs

# num <- round(10^(seq(2,log10(nrow(R)),length=12)))
num <- round(10^(seq(2,4.8,0.2)))
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(R), n)))
# sample_index <- sapply(1:100, function(i) sample(nrow(R), n))
# kdv_ada_lasso_pareto <- function(sam, ...){
#   coeff_full <- ada_lasso_pareto_f(R[sam,], Ut[sam,], ...)[[1]]
#   names(coeff_full) <- rhs
#   coeff <- coeff_full[which(coeff_full!=0)]
#   return(coeff)
# }

kdv_rudy <- function(sam, lam=10^-5, d_tol=5, ...){
  # rudy
  w = TrainSTRidge(R[sam,], as.matrix(Ut)[sam,], lam, d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- rhs[which(w!=0)]
  return(coef.rudy)
}

d_tols <- 1:10
d_tols <- seq(10,100,10)
d_tols <- c(0.2, 2, 10)
system.time(
  outs_rudy <- mclapply(d_tols, function(tol){
    outs <- lapply(1:length(num), function(i){
      rudy_out_list <- list()
      for(j in 1:100){
        rudy_out_list[[j]] <- kdv_rudy(sample_index[[i]][[j]], d_tol=tol)
      }
      return(rudy_out_list)
    })
    names(outs) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
    return(outs)
  }, mc.cores=10)
)


dens_fun_kdv <- function(noise_coeff, true_terms){
  dens <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(dens)
}
kdv_true_terms <- c("uu_{x}", "u_{xxx}")
success_rate_rudy_list <- lapply(outs_rudy, function(k) sapply(seq_along(num), function(i) dens_fun_kdv(k[[i]], kdv_true_terms)))

save(outs_rudy, success_rate_rudy_list,
     file=sprintf('Tests/Outputs/kdv_SG_noiseless_seed_%s_success_rate_rudy_d_thred.RData', seed))

rm(list=ls())
