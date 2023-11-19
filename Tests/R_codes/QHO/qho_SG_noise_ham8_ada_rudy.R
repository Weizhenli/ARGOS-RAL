setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate);library(doParallel)
library(R.matlab);library(SMFilter);library(signal)
cv <- import('cv2')
seed = 10
set.seed(seed)
source('Functions/ada_lasso_pareto.R')

# qho ----------------------
source('Functions/all_functions.R')

source("pde_solver_data/qho_solver.R")
qho.x = xx
n <- length(xx)
m <- length(tt)

build_qho_Theta <- function(u, dx, dt, sg_para=NULL, deg_min=3, deg_max=6, width_max=10, noise=0, diff_o=3, poly_o=3,
                            ref_data='Gblur', ...){ 
  u <- as.matrix(u)
  if(noise > 0){
    un <- u + rnorm(nrow(u)*ncol(u),sd=0.01/sqrt(2)*sd(Re(u))) + 1i*rnorm(nrow(u)*ncol(u),sd=0.01/sqrt(2)*sd(Im(u)))
  }else{
    un <- u
    #un_gblur <- u
  }
  if(ref_data == 'Gblur'){
    un_gblur <- cv$GaussianBlur(Re(un),as.integer(c(3,3)),0) + 1i*cv$GaussianBlur(Im(un),as.integer(c(3,3)),0)
  }else{
    un_gblur <- un
  }
  
  n=nrow(un);m=ncol(un)
  if(is.null(sg_para)){
    deg=deg_min:deg_max
    mse_t = mse_x <- matrix(Inf,nrow=deg_max,ncol=width_max)
    for(d in deg){
      width_min <- (d+3-(d %% 2)-1)/2
      for(w in width_min:width_max){
        un_sg <- matrix(0,nrow=n,ncol=m)
        for(i in 1:n){
          un_sg[i,]<-abs(sgolayfilt(Re(un[i,]),p=d,n=2*w+1)+1i*sgolayfilt(Im(un[i,]),p=d,n=2*w+1))
        }
        mse_t[d,w] <- sum((un_sg-abs(un_gblur))^2)/(nrow(un_sg)*ncol(un_sg))
        un_sg <- matrix(0,nrow=n,ncol=m)
        for(i in 1:m){
          un_sg[,i]<-abs(sgolayfilt(Re(un[,i]),p=d,n=2*w+1)+1i*sgolayfilt(Im(un[,i]),p=d,n=2*w+1))
        }
        mse_x[d,w] <- sum((un_sg-abs(un_gblur))^2)/(nrow(un_sg)*ncol(un_sg))
      }
    }
    deg_x <- which.min(apply(mse_x,1,min)) # index of row
    width_x <- which.min(apply(mse_x,2,min)) # index of column
    
    deg_t <- which.min(apply(mse_t,1,min)) # index of row
    width_t <- which.min(apply(mse_t,2,min)) # index of column
  }else{
    deg_x <- sg_para$deg_x # index of row
    width_x <- sg_para$width_x # index of column
    deg_t <- sg_para$deg_t # index of row
    width_t <- sg_para$width_t # index of column
  }
  
  n2 = n-2*width_t
  m2 = m-2*width_x
  potential <- matrix(rep(0.5*qho.x[width_x:(m-width_x-1)]^2,each=n2),nrow=n2)
  dim(potential) <- c(m2*n2, 1)
  
  utn <- matrix(complex(), nrow=n2, ncol=m2)
  uxn <- matrix(complex(), nrow=n2, ncol=m2)
  uxxn <- matrix(complex(), nrow=n2, ncol=m2)
  uxxxn <- matrix(complex(), nrow=n2, ncol=m2)
  
  x_deriv <- list()
  for(j in 1:diff_o){
    x_deriv <- c(x_deriv, list(matrix(NA, nrow=n2, ncol=m2)))
  }
  
  for(i in 1:m2){
    utn[,i] = PolyDiff(Re(un_gblur[,i+width_x]), dx=dt, deg=deg_t, width=width_t, diff=1)
    utn[,i] = utn[,i] + 1i*PolyDiff(Im(un_gblur[,i+width_x]), dx=dt, deg=deg_t, width=width_t, diff=1)
  }
  for(i in 1:n2){
    x_derivatives = PolyDiff(Re(un_gblur[i+width_t,]), dx, deg=deg_x, diff=diff_o, width=width_x)
    x_derivatives = x_derivatives + 1i*PolyDiff(Im(un_gblur[i+width_t,]), dx, deg=deg_x, diff=diff_o, width=width_x)
    for(j in 1:diff_o){
      x_deriv[[j]][i,] <- x_derivatives[,j]
    }
  }
  
  dim(utn) <- c(n2*m2, 1)
  for(j in 1:diff_o){
    dim(x_deriv[[j]]) <- c(n2*m2, 1)
  }
  
  X_ders <- cbind(rep(1,n2*m2), as.matrix(as.data.frame(x_deriv, col.names = 1:diff_o)))
  un_gblur.copy <- matrix(un_gblur[(width_t+1):(n-width_t),(width_x+1):(m-width_x)],nrow=n2*m2,ncol=1)
  un_gblur.abs <- matrix(abs(un_gblur[(width_t+1):(n-width_t),(width_x+1):(m-width_x)]),nrow=n2*m2,ncol=1)
  X_data <- cbind(un_gblur.copy,un_gblur.abs,potential)
  derivatives_description <- c('','u_{x}','u_{xx}', 'u_{xxx}')
  data_description <- c('u','|u|','V')
  X1n <- build_Theta(X_data, X_ders, derivatives_description, poly_o, data_description)
  qho_reg_data <- transfer(utn,X1n) # transfer complex data to regressionable
  Ut <- qho_reg_data[, 1]
  R <- qho_reg_data[, -1]
  rhs <- colnames(qho_reg_data)[-1]
  
  return(list(ut = Ut, Theta = R, rhs = rhs))
}

snr_db_seq <- c(seq(40, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(eta_i){
  lapply(1:100, function(x){
    rnorm(m*n, sd=eta_i/sqrt(2)*sd(Re(qho.usol))) + 1i*rnorm(m*n, sd=eta_i/sqrt(2)*sd(Im(qho.usol)))
  })
})

qho_all_density_noise_fun <- function(eta_i, u, num=100, lam=10^-5, multi_cores=NULL, ...){
  noise_list1 <- list(rudy_0.2=list(), rudy_2=list(), rudy_5=list(), rudy_10=list(), ada_lasso=list())
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  multi_cores <- ifelse(is.null(multi_cores),1,multi_cores)
  noise_list0 <- mclapply(1:num, function(i){
    noise_list <- list(rudy_0.2=list(), rudy_2=list(), rudy_5=list(), rudy_10=list(), ada_lasso=list())
    # qho_un <- u + rnorm(m*n, sd=eta/sqrt(2)*sd(Re(u))) + 1i*rnorm(m*n, sd=eta/sqrt(2)*sd(Im(u))) 
    qho_un <- u + u_noise[[eta_i]][[i]]
    candidate_library <- build_qho_Theta(qho_un, dx=dx, dt=dt,noise=0)
    Ut <- candidate_library$ut
    R <- candidate_library$Theta
    rhs <- candidate_library$rhs
    colnames(R) <- rhs
    # rudy
    w = TrainSTRidge(R, Ut, lam, d_tol=0.2)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    noise_list$rudy_0.2 <- coef.rudy
    
    w = TrainSTRidge(R, Ut, lam, d_tol=2)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    noise_list$rudy_2 <- coef.rudy
    
    w = TrainSTRidge(R, Ut, lam, d_tol=5)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    noise_list$rudy_5 <- coef.rudy
    
    w = TrainSTRidge(R, Ut, lam, d_tol=10)
    coef.rudy <- w[which(w!=0)]
    names(coef.rudy) <- rhs[which(w!=0)]
    noise_list$rudy_10 <- coef.rudy
    
    # reweighted ada lasso 
    coeff_full_ada <- ada_lasso_pareto_f_complex(R, Ut, ...)[[1]]
    names(coeff_full_ada) <- rhs
    coeff_ada <- coeff_full_ada[which(coeff_full_ada!=0)]
    noise_list$ada_lasso <- coeff_ada
    return(noise_list)
  }, mc.cores=multi_cores)
  noise_list1$rudy_0.2 <- lapply(noise_list0, function(k) k$rudy_0.2)
  noise_list1$rudy_2 <- lapply(noise_list0, function(k) k$rudy_2)
  noise_list1$rudy_5 <- lapply(noise_list0, function(k) k$rudy_5)
  noise_list1$rudy_10 <- lapply(noise_list0, function(k) k$rudy_10)
  noise_list1$ada_lasso <- lapply(noise_list0, function(k) k$ada_lasso)
  return(noise_list1)
}

sample_size <- 100
library(doParallel) # start parallel computation
no_cores <- 20 # number of cpu cores that will be used 

system.time(
  # qho_noise_all_results <- mclapply(eta, qho_all_density_noise_fun, u=qho.usol, num=sample_size, mc.cores=no_cores)
  qho_noise_all_results <- lapply(seq_along(eta), qho_all_density_noise_fun, u=qho.usol, num=sample_size, multi_cores=no_cores)
)
save(qho_noise_all_results, file=sprintf('/nobackup/cfzh32/for_r/SG/qho/qho_SG_noise_seed_%s_samp_%s_snr_ada_rudy_back.RData', seed, sample_size))

qho_ada_lasso_noise <- lapply(qho_noise_all_results, function(x) x$ada_lasso)
qho_rudy_0.2_noise <- lapply(qho_noise_all_results, function(x) x$rudy_0.2)
qho_rudy_2_noise <- lapply(qho_noise_all_results, function(x) x$rudy_2)
qho_rudy_5_noise <- lapply(qho_noise_all_results, function(x) x$rudy_5)
qho_rudy_10_noise <- lapply(qho_noise_all_results, function(x) x$rudy_10)



qho_true_term <- c("iu_{xx}", "iuV")
qho_bernoulli_fun <- function(qho_data){
  # qho_noise <- qho_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(qho_data), function(i){
    match_index <- match(names(qho_data[[i]]), qho_true_term)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

qho_density_noise_rudy_0.2 <- sapply(qho_rudy_0.2_noise, qho_bernoulli_fun)
qho_density_noise_rudy_2 <- sapply(qho_rudy_2_noise, qho_bernoulli_fun)
qho_density_noise_rudy_5 <- sapply(qho_rudy_5_noise, qho_bernoulli_fun)
qho_density_noise_rudy_10 <- sapply(qho_rudy_10_noise, qho_bernoulli_fun)
qho_density_noise_ada_AIC <- sapply(qho_ada_lasso_noise, qho_bernoulli_fun)

save(qho_ada_lasso_noise, qho_noise_all_results,
     qho_density_noise_back, qho_density_noise_ada_AIC,
     qho_rudy_0.2_noise, qho_density_noise_rudy_0.2,
     qho_rudy_2_noise, qho_density_noise_rudy_2,
     qho_rudy_5_noise, qho_density_noise_rudy_5,
     qho_rudy_10_noise, qho_density_noise_rudy_10,
     file=sprintf('Tests/Outputs/qho_SG_noise_seed_%s_samp_%s_snr_ada_rudy_back.RData', seed, sample_size))
rm(list=ls())


