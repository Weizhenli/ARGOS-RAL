setwd('/mnt/d/GitHub/ARGOS-RAL/')
rm(list=ls())
library(glmnet);library(leaps);library(reticulate) 
library(R.matlab);library(SMFilter);library(signal)
library(doParallel)
cv <- import('cv2')
source('Functions/all_functions.R')
source('Functions/ada_lasso_pareto.R')
source('Functions/SG_opt_function.R')

diff_o <- as.numeric(Sys.getenv("d"))
poly_o <- as.numeric(Sys.getenv("p"))

out_list = F_test_list <- list()
set.seed(1)
dt=0.01; dx=0.02
n=2000000; nrow=2000; mean=0;
sd=as.numeric(Sys.getenv("sd"))
## Use for loop ---------------------------------
# for(i in 1:100){
#   test_dataset <- matrix(rnorm(n, mean, sd), nrow=nrow)
#   test_dataset <- as.data.frame(test_dataset)
#   # print(test_dataset[1:5,1:5])
#   # image(as.matrix(test_dataset[,-1]))
#   candidate_library <- ASG_build_library(test_dataset,dx,dt,noise=0, diff_o=diff_o, poly_o=poly_o)
#   candidate_library$rhs[1] <- '1'
#   sample_index <- sample(1:nrow(candidate_library$Theta),1e+4)
#   # candidate_library$Theta[,1] <- rnorm(nrow(candidate_library$Theta),1,1e+15)
#   test_random_data <- cbind.data.frame(y=candidate_library$ut, X=candidate_library$Theta)[sample_index,]
#   colnames(test_random_data) <- c('ut',candidate_library$rhs)
#   coeff_random <- ada_lasso_pareto_f(test_random_data[,-1], test_random_data[,1]) # RAL
#   ARGOS_active <- which(coeff_random[[1]]!=0)
#   ARGOS_lm_model <- lm(ut~0+., data = test_random_data[,c(0,ARGOS_active)+1])
#   ARGOS_lm <- summary(ARGOS_lm_model)
#   ARGOS_coef <- coeff_random[[2]]
#   names(ARGOS_coef) <- as.character(ARGOS_active)
#   w_all_random = TrainSTRidge(test_random_data[,-1], test_random_data[,1], 1e-5, 2) # TSTRidge
#   TSTRidge_active <- which(w_all_random!=0)
#   TSTRidge_model <- lm(ut~0+., data = test_random_data[,c(0,TSTRidge_active)+1])
#   TSTRidge_lm <- summary(TSTRidge_model)
#   w_random <- w_all_random[w_all_random!=0]
#   names(w_random) <- as.character(TSTRidge_active)
#   F_test_list[[i]] <- list(pf(ARGOS_lm$fstatistic[1], ARGOS_lm$fstatistic[2], ARGOS_lm$fstatistic[3], lower.tail = FALSE),
#                            pf(TSTRidge_lm$fstatistic[1], TSTRidge_lm$fstatistic[2], TSTRidge_lm$fstatistic[3], lower.tail = FALSE))
#   
#   out_list[[i]] <- list(ARGOS=ARGOS_coef, PDE_FIND=w_random)
#   print(list(ARGOS=coeff_random[[1]][coeff_random[[1]]!=0], PDE_FIND=w_random))
# }

## parallel running ----------------------------------
out_list_all <- mclapply(1:100, function(i){
  test_dataset <- matrix(rnorm(n, mean, sd), nrow=nrow)
  test_dataset <- as.data.frame(test_dataset)
  candidate_library <- ASG_build_library(test_dataset,dx,dt,noise=0, diff_o=diff_o, poly_o=poly_o)
  candidate_library$rhs[1] <- '1'
  sample_index <- sample(1:nrow(candidate_library$Theta),1e+4)
  test_random_data <- cbind.data.frame(y=candidate_library$ut, X=candidate_library$Theta)[sample_index,]
  colnames(test_random_data) <- c('ut',candidate_library$rhs)
  coeff_random <- ada_lasso_pareto_f(test_random_data[,-1], test_random_data[,1]) # RAL
  ARGOS_active <- which(coeff_random[[1]]!=0)
  ARGOS_lm_model <- lm(ut~0+., data = test_random_data[,c(0,ARGOS_active)+1])
  ARGOS_lm <- summary(ARGOS_lm_model)
  ARGOS_coef <- coeff_random[[2]]
  names(ARGOS_coef) <- as.character(ARGOS_active)
  w_all_random = TrainSTRidge(test_random_data[,-1], test_random_data[,1], 1e-5, 2) # TSTRidge
  TSTRidge_active <- which(w_all_random!=0)
  TSTRidge_model <- lm(ut~0+., data = test_random_data[,c(0,TSTRidge_active)+1])
  TSTRidge_lm <- summary(TSTRidge_model)
  w_random <- w_all_random[w_all_random!=0]
  names(w_random) <- as.character(TSTRidge_active)
  F_test <- list(pf(ARGOS_lm$fstatistic[1], ARGOS_lm$fstatistic[2], ARGOS_lm$fstatistic[3], lower.tail = FALSE),
                           pf(TSTRidge_lm$fstatistic[1], TSTRidge_lm$fstatistic[2], TSTRidge_lm$fstatistic[3], lower.tail = FALSE))
  
  out <- list(ARGOS=ARGOS_coef, PDE_FIND=w_random)
  print(list(ARGOS=coeff_random[[1]][coeff_random[[1]]!=0], PDE_FIND=w_random))
  return(list(out = out, F_test = F_test))
}, mc.cores=20)

out_list <- lapply(out_list_all, function(x) x[[1]])
F_test_list <- lapply(out_list_all, function(x) x[[2]])

ARGOS <- lapply(out_list, function(x) x[[1]])
PDE_FIND <- lapply(out_list, function(x) x[[2]])
ARGOS_fp <- lapply(F_test_list, function(x) x[[1]])
PDE_FIND_fp <- lapply(F_test_list, function(x) x[[2]])
terms_names <- candidate_library$rhs

save(ARGOS,ARGOS_fp, PDE_FIND,PDE_FIND_fp,terms_names, file = sprintf('Tests/Outputs/white_noise_test_out_p%s_d%s_sd%s.RData',poly_o,diff_o,sd))
