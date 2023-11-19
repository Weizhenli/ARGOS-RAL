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
    #print(c(deg_x,width_x,deg_t,width_t))
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

candidate_library <- build_qho_Theta(qho.usol, dx=dx, dt=dt, noise=0)
Ut <- candidate_library$ut
R <- candidate_library$Theta
rhs <- candidate_library$rhs
colnames(R) <- rhs

num <- 10^(seq(2.2,5.2,0.2))
sample_index <- lapply(num, function(n) lapply(1:100, function(i){
  index <- sample(nrow(R)/2, n)
  return(c(index*2-1, index*2))
}))


qho_rudy_density_noiseless_fun <- function(sam, lam=10^-5, d_tol=10, ...){
  # rudy
  w = TrainSTRidge(R[sam,], as.matrix(Ut)[sam,], lam, d_tol)
  coef.rudy <- w[which(w!=0)]
  names(coef.rudy) <- rhs[which(w!=0)]
  return(coef.rudy)
}

qho_ada_density_noiseless_fun <- function(sam, ...){
  # reweighted ada lasso 
  coeff_full_ada <- ada_lasso_pareto_f_complex(R[sam,], as.matrix(Ut)[sam,], ...)[[1]]
  names(coeff_full_ada) <- rhs
  coeff_ada <- coeff_full_ada[which(coeff_full_ada!=0)]
  return(coeff_ada)
}

library(doParallel) # start parallel computation
no_cores <- 25 # number of cpu cores that will be used 

qho_noiseless_all_results <- list()
system.time(
  for(i in 1:length(num)){
    noiseless_list <- mclapply(1:100, function(j){
      rudy_0.2 <- qho_rudy_density_noiseless_fun(sample_index[[i]][[j]],d_tol=0.2)
      rudy_2 <- qho_rudy_density_noiseless_fun(sample_index[[i]][[j]],d_tol=2)
      rudy_5 <- qho_rudy_density_noiseless_fun(sample_index[[i]][[j]],d_tol=5)
      rudy_10 <- qho_rudy_density_noiseless_fun(sample_index[[i]][[j]],d_tol=10)
      ada_lasso <- qho_ada_density_noiseless_fun(sample_index[[i]][[j]])
      return(list(rudy_0.2 = rudy_0.2, rudy_2 = rudy_2, rudy_5 = rudy_5, rudy_10 = rudy_10, ada_lasso=ada_lasso))
    }, mc.cores=no_cores)
    qho_noiseless_all_results[[i]] <- noiseless_list
    save(qho_noiseless_all_results,
         file=sprintf('Tests/Outputs/qho_SG_noiseless_seed_%s_n_ada_rudy_back.RData', seed))
    print(i)
  }
)

# save(qho_noiseless_all_results, file=sprintf('Tests/Outputs/qho_SG_noiseless_seed_%s_n_ada_rudy_back.RData', seed))

qho_ada_lasso_noiseless <- lapply(qho_noiseless_all_results, function(x) lapply(x, function(y) y$ada_lasso))
qho_rudy_0.2_noiseless <- lapply(qho_noiseless_all_results, function(x) lapply(x, function(y) y$rudy_0.2))
qho_rudy_2_noiseless <- lapply(qho_noiseless_all_results, function(x) lapply(x, function(y) y$rudy_2))
qho_rudy_5_noiseless <- lapply(qho_noiseless_all_results, function(x) lapply(x, function(y) y$rudy_5))
qho_rudy_10_noiseless <- lapply(qho_noiseless_all_results, function(x) lapply(x, function(y) y$rudy_10))

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

qho_density_noiseless_rudy_0.2 <- sapply(qho_rudy_0.2_noiseless, qho_bernoulli_fun)
qho_density_noiseless_rudy_2 <- sapply(qho_rudy_2_noiseless, qho_bernoulli_fun)
qho_density_noiseless_rudy_5 <- sapply(qho_rudy_5_noiseless, qho_bernoulli_fun)
qho_density_noiseless_rudy_10 <- sapply(qho_rudy_10_noiseless, qho_bernoulli_fun)
qho_density_noiseless_ada_AIC <- sapply(qho_ada_lasso_noiseless, qho_bernoulli_fun)

save(qho_ada_lasso_noiseless, qho_noiseless_all_results,
     qho_density_noiseless_back, qho_density_noiseless_ada_AIC,
     qho_rudy_0.2_noiseless,qho_density_noiseless_rudy_0.2,
     qho_rudy_2_noiseless,qho_density_noiseless_rudy_2,
     qho_rudy_5_noiseless,qho_density_noiseless_rudy_5,
     qho_rudy_10_noiseless,qho_density_noiseless_rudy_10,
     file=sprintf('Tests/Outputs/qho_SG_noiseless_seed_%s_n_ada_rudy_back.RData', seed))
rm(list=ls())


