setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
seed = 100
set.seed(seed)
crit = 'AIC'

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

n <- length(burgers.t) # 101
m <- length(burgers.x) # 256

bur_true_terms <- c("uu_{x}", "u_{xx}")

# create noise 
snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)

u_noise <- lapply(eta, function(noise){
  lapply(1:100, function(x){
    matrix(rnorm(m*n, 0, noise*sd(burgers.usol)), nrow=m, ncol=n)
  })
})

density_noise_fun <- function(eta, u, num=100, seed=100, ...){
  # set.seed(seed)
  noise_list <- list()
  u <- as.matrix(u)
  m <- nrow(u)
  n <- ncol(u)
  for(i in 1:num){
    bur_un <- u + matrix(rnorm(m*n, 0, eta*sd(u)), nrow=m, ncol=n) # in linux it has problem
    candidate_library <- ASG_build_library(bur_un,dx,dt,noise=0)
    Ut <- candidate_library$ut
    R <- candidate_library$Theta
    rhs <- candidate_library$rhs
    coeff_full <- ada_lasso_pareto_f(R, Ut, ...)[[1]]
    names(coeff_full) <- rhs
    coeff <- coeff_full[which(coeff_full!=0)]
    noise_list[[i]] <- coeff
  }
  return(noise_list)
}

sample_size <- 100
library(doParallel) # start parallel computation
no_cores <- 22 # number of cpu cores that will be used 
cl <- makeCluster(no_cores, type="FORK")
# system.time(bur_all_noise <- parLapply(cl, 0:10*0.01, density_noise_fun, u=burgers.usol, num=sample_size, seed=seed))
system.time(bur_ada_lasso_pareto_noise <- parLapply(cl, eta, density_noise_fun, u=burgers.usol, num=sample_size, seed=seed, crit=crit))
stopCluster(cl)
save(bur_ada_lasso_pareto_noise, file=sprintf('Tests/Outputs/bur_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto_%s.RData', seed, sample_size,crit))

bur_density_fun <- function(i, bur_noise_data){
  bur_noise <- bur_noise_data[[i]]
  bernoulli_noise <- sapply(seq_along(bur_noise), function(i){
    match_index <- match(names(bur_noise[[i]]), bur_true_terms)
    if(length(match_index) == 2 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(bernoulli_noise)
}

bur_density_noise_ada_lasso_pareto <- lapply(list(bur_ada_lasso_pareto_noise), function(q) sapply(seq_along(q), bur_density_fun, bur_noise_data=q))

save(bur_ada_lasso_pareto_noise, bur_density_noise_ada_lasso_pareto,
     file=sprintf('Tests/Outputs/bur_SG_noise_seed_%s_samp_%s_snr_ada_lasso_pareto_%s.RData', seed, sample_size,crit))
rm(list=ls())
