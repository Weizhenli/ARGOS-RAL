setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)

source_python('pde_solver_data/ode_int.py')
source('Functions/all_functions.R')
cv <- import('cv2')
source('Functions/ada_lasso_pareto.R')

source_python('Tests/Outputs/more_pde_test/Keller_Segel_model.py')

u <- t(u_u)
v <- t(v_u)

# Only finite difference works to identify the equation
poly_o=3
## finite difference
ut0 <- t(apply(u, 1, FiniteDiff, dx=dt, d=1))
ux <- apply(u, 2, FiniteDiff, dx=dx, d=1)
uxx <- apply(u, 2, FiniteDiff, dx=dx, d=2)
uxxx <- apply(u, 2, FiniteDiff, dx=dx, d=3)
vt0 <- t(apply(v, 1, FiniteDiff, dx=dt, d=1))
vx <- apply(v, 2, FiniteDiff, dx=dx, d=1)
vxx <- apply(v, 2, FiniteDiff, dx=dx, d=2)
vxxx <- apply(v, 2, FiniteDiff, dx=dx, d=3)

ut <- as.matrix(as.vector(ut0))
vt <- as.matrix(as.vector(vt0))
u_ders <- cbind(1, as.vector(ux), as.vector(uxx), as.vector(uxxx))
v_ders <- cbind(1, as.vector(vx), as.vector(vxx), as.vector(vxxx))
X_data <- cbind(as.vector(u),as.vector(v))

## candidate lib
derivatives_description <- c('u_{x}','u_{xx}', 'u_{xxx}')
data_description <- c('u','v')
X1n <- build_Theta(X_data, v_ders, c('','v_{x}','v_{xx}', 'v_{xxx}'), poly_o, data_description)
X1n <- X1n[,-1]
rhs <- names(X1n)
# rhs[1] <- '1'
X1n2 <- build_Theta(X1n, u_ders, c('','u_{x}','u_{xx}', 'u_{xxx}'), 1, rhs)
rhs2 <- names(X1n2)


num <- 10^(seq(2.4,5,0.2)) # minimum n should be 160
sample_index <- lapply(num, function(n) lapply(1:100, function(i) sample(nrow(X1n2), n)))
Keller_ada_lasso_pareto <- function(sam, ...){
  coeff_full_ut <- ada_lasso_pareto_f(X1n2[sam,], ut[sam,], ...)[[1]]
  names(coeff_full_ut) <- rhs2
  coeff_ut <- coeff_full_ut[which(coeff_full_ut!=0)]
  coeff_full_vt <- ada_lasso_pareto_f(X1n2[sam,], vt[sam,], ...)[[1]]
  names(coeff_full_vt) <- rhs2
  coeff_vt <- coeff_full_vt[which(coeff_full_vt!=0)]
  return(list(coeff_ut,coeff_vt))
}

Keller_stepwise_backward <- function(sam, ...){
  out_ut <- path_stepwise(X1n2[sam,], ut[sam,], 3, ...)[[1]][[1]]
  out_vt <- path_stepwise(X1n2[sam,], vt[sam,], 3, ...)[[1]][[1]]
  return(list(out_ut,out_vt))
}

Keller_rudy <- function(sam, lam=10^-5, d_tol=1, ...){
  w_ut = TrainSTRidge(X1n2[sam,], ut[sam,], lam=lam, d_tol=d_tol)
  coef.rudy_ut <- w_ut[which(w_ut!=0)]
  names(coef.rudy_ut) <- rhs2[which(w_ut!=0)]
  w_vt = TrainSTRidge(X1n2[sam,], vt[sam,], lam=lam, d_tol=d_tol)
  coef.rudy_vt <- w_vt[which(w_vt!=0)]
  names(coef.rudy_vt) <- rhs2[which(w_vt!=0)]
  return(list(coef.rudy_ut,coef.rudy_vt))
}

d_tols <- c(0.2, 2, 10)
library(doParallel)
# rudy 
system.time(
  outs_rudy_list <- mclapply(d_tols, function(tol){
    outs <- mclapply(1:length(num), function(i){
      rudy_out_list <- list()
      for(j in 1:100){
        rudy_out_list[[j]] <- Keller_rudy(sort(sample_index[[i]][[j]]), d_tol=tol)
      }
      return(rudy_out_list)
    }, mc.cores=10)
    names(outs) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
    return(outs)
  }, mc.cores=3)
)

# ada-lasso
system.time(outs_ada_lasso_pareto_success_rate <- mclapply(1:length(num), function(i){
  ada_lasso_pareto_out_list <- list()
  ada_lasso_pareto_out_list <- mclapply(1:100, function(j){
    Keller_ada_lasso_pareto(sort(sample_index[[i]][[j]]))
  }, mc.cores=50)
  return(ada_lasso_pareto_out_list)
}, mc.cores=1))
names(outs_ada_lasso_pareto_success_rate) <- sapply(1:length(num), function(j) paste0('n=',round(j)))
# back
system.time(outs_back_success_rate <- mclapply(1:length(num), function(i){
  back_out_list <- list()
  for(j in 1:100){
    back_out_list[[j]] <- Keller_stepwise_backward(sort(sample_index[[i]][[j]]))
  }
  return(back_out_list)
}, mc.cores=length(num)))
names(outs_back_success_rate) <- sapply(1:length(num), function(j) paste0('n=',round(j)))

save(outs_rudy_list, outs_ada_lasso_pareto_success_rate, outs_back_success_rate, 
     file = 'Tests/Outputs/keller_segel_FD_noiseless_ham8_ada_rudy_back.RData')

dens_fun_KS <- function(noise_coeff, true_terms, which){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index_ut <- match(names(noise_coeff[[i]][[1]]), true_terms[[1]])
    match_index_vt <- match(names(noise_coeff[[i]][[2]]), true_terms[[2]])
    match_TF_all <- length(match_index_ut) == 3 & all(!is.na(match_index_ut)) & 
      length(match_index_vt) == 3 & all(!is.na(match_index_vt))
    match_TF_ut <- length(match_index_ut) == 3 & all(!is.na(match_index_ut))
    match_TF_vt <- length(match_index_vt) == 3 & all(!is.na(match_index_vt))
    if(which=="ut"){
      if(match_TF_ut){
        return(1)
      }else{
        return(0)
      }
    }else if(which=="vt"){
      if(match_TF_vt){
        return(1)
      }else{
        return(0)
      }
    }else{
      if(match_TF_all){
        return(1)
      }else{
        return(0)
      }
    }
  })
  return(dens_noise)
}
true_term_KS <- list(u_t = c("u_{xx}", "v_{x}u_{x}", "uv_{xx}"),
                     v_t = c("v_{xx}", "u", "v"))


dens_our <- sapply(seq_along(num), function(i) dens_fun_KS(outs_ada_lasso_pareto_success_rate[[i]], true_term_KS, "all"))
dens_rudy_list <- lapply(outs_rudy_list, function(k) sapply(seq_along(num), function(i) dens_fun_KS(k[[i]], true_term_KS, "all"))) 
dens_back <- sapply(seq_along(num), function(i) dens_fun_KS(outs_back_success_rate[[i]], true_term_KS, "all"))


save(outs_rudy_list, outs_ada_lasso_pareto_success_rate, outs_back_success_rate, 
     dens_our, dens_rudy_list, dens_back, 
     file = 'Tests/Outputs/keller_segel_FD_noiseless_ham8_ada_rudy_back.RData')
