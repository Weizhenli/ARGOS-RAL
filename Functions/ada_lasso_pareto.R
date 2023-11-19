library(glmnet)
source('/nobackup/cfzh32/for_r/functions/opt_point_L-curve_corner.R')
library(reticulate)
cv <- import('cv2')


################################################################
## transform the algorithms oncultrera_simple_2020
## A simple algorithm to find the L-curve corner in the regularization of ill-posed inverse problems
## to R codes

l_curve_P0 <- function(glmnet_obj, y, x, weights){
  # glmnet_obj is a glmnet object
  residuals = as.numeric(y)- as.matrix(x) %*% glmnet_obj$beta
  xi <- log(apply(residuals,2,norm, type='2'))
  eta <- log(as.numeric(weights %*% abs(glmnet_obj$beta)))
  return(cbind(xi,eta))
}


l_curve_P <- function(glmnet_obj, y, x, weights, lambda){
  # glmnet_obj is a glmnet object
  glmnet_coeff <- coef(glmnet_obj, s = lambda)[-1]
  residuals = as.numeric(y)- as.matrix(x) %*% glmnet_coeff
  xi <- log(apply(residuals,2,norm, type='2'))
  alpha = glmnet_obj$call$alpha
  eta <- log(as.numeric(weights %*% (alpha*abs(glmnet_coeff)+(1-alpha)/2*glmnet_coeff^2)))
  return(c(xi,eta))
}


menger <- function(pj, pk, pl){
  pjpk <- (pk[1] - pj[1])^2 + (pk[2] - pj[2])^2
  pkpl <- (pl[1] - pl[1])^2 + (pl[2] - pk[2])^2
  pkpl <- (pj[1] - pl[1])^2 + (pj[2] - pl[2])^2
  C_k <- 2*(pj[1]*pk[2]+pk[1]*pl[2]+pl[1]*pj[2]-pj[1]*pl[2]-pk[1]*pj[2]-pl[1]*pk[2])/sqrt(pjpk*pkpl*pkpl)
  return(C_k)
}

l_curve_corner_search <- function(glmnet_obj, y, x, weights, lambda1,lambda4, stop_threshold = 1e-7){
  phi <- (1 + sqrt(5)) / 2
  x1 <- log10(lambda1) 
  x4 <- log10(lambda4)
  x2 <- (x4 + phi * x1) / (1+phi) 
  x3 <- x1 + x4 - x2
  lambda2 <- 10^(x2)
  lambda3 <- 10^(x3)
  l <- c(lambda1, lambda2, lambda3, lambda4)
  p <- list()
  for(i in 1:4){
    p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
  }
  if(all(p[[1]] == p[[4]])){
    stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
  }
  
  while(all(p[[1]] == p[[3]])){
    lambda1 <- lambda3
    x1 <- log10(lambda1)
    x4 <- log10(lambda4)
    x2 <- (x4 + phi * x1) / (1+phi) 
    x3 <- x1 + x4 - x2
    lambda2 <- 10^(x2)
    lambda3 <- 10^(x3)
    l <- c(lambda1, lambda2, lambda3, lambda4)
    p <- list()
    for(i in 1:4){
      p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
    }
    if(all(p[[1]] == p[[4]])){
      stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
    }
  }
  while(all(p[[2]] == p[[4]])){
    lambda4 <- lambda2
    x1 <- log10(lambda1)
    x4 <- log10(lambda4)
    x2 <- (x4 + phi * x1) / (1+phi) 
    x3 <- x1 + x4 - x2
    lambda2 <- 10^(x2)
    lambda3 <- 10^(x3)
    l <- c(lambda1, lambda2, lambda3, lambda4)
    p <- list()
    for(i in 1:4){
      p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
    }
    if(all(p[[1]] == p[[4]])){
      stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
    }
  }
  while(all(p[[1]] == p[[2]])){
    lambda1 <- lambda2
    x1 <- log10(lambda1)
    x4 <- log10(lambda4)
    x2 <- (x4 + phi * x1) / (1+phi) 
    x3 <- x1 + x4 - x2
    lambda2 <- 10^(x2)
    lambda3 <- 10^(x3)
    l <- c(lambda1, lambda2, lambda3, lambda4)
    p <- list()
    for(i in 1:4){
      p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
    }
    if(all(p[[1]] == p[[4]])){
      stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
    }
  }
  while(all(p[[2]] == p[[3]])){
    lambda1 <- lambda2
    x1 <- log10(lambda1)
    x4 <- log10(lambda4)
    x2 <- (x4 + phi * x1) / (1+phi) 
    x3 <- x1 + x4 - x2
    lambda2 <- 10^(x2)
    lambda3 <- 10^(x3)
    l <- c(lambda1, lambda2, lambda3, lambda4)
    p <- list()
    for(i in 1:4){
      p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
    }
    if(all(p[[1]] == p[[4]])){
      stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
    }
  }
  while(all(p[[3]] == p[[4]])){
    lambda4 <- lambda3
    x1 <- log10(lambda1)
    x4 <- log10(lambda4)
    x2 <- (x4 + phi * x1) / (1+phi) 
    x3 <- x1 + x4 - x2
    lambda2 <- 10^(x2)
    lambda3 <- 10^(x3)
    l <- c(lambda1, lambda2, lambda3, lambda4)
    p <- list()
    for(i in 1:4){
      p[[i]] <- l_curve_P(glmnet_obj, y, x, weights, l[i])
    }
    if(all(p[[1]] == p[[4]])){
      stop(paste('The coefficients from',l[1],'and',l[4],'are the same.'))
    }
  }
  
  while((lambda4-lambda1)/lambda4 >= stop_threshold){
    # print((lambda4-lambda1)/lambda4)
    C2 <- menger(p[[1]],p[[2]],p[[3]])
    C3 <- menger(p[[2]],p[[3]],p[[4]])
    while(C3 <= 0){
      lambda4 <- lambda3; p[[4]] <- p[[3]]
      lambda3 <- lambda2; p[[3]] <- p[[2]]
      x1 <- log10(lambda1) 
      x4 <- log10(lambda4)
      x2 <- (x4 + phi * x1) / (1+phi) 
      x3 <- x1 + x4 - x2
      lambda2 <- 10^(x2)
      p[[2]] <- l_curve_P(glmnet_obj, y, x, weights, lambda2)
      C3 <- menger(p[[2]],p[[3]],p[[4]])
      # C3 may not converge = -Inf
      if(abs(C3) == Inf){
        break
      }
    }
    if(C2 > C3){
      lambda <- lambda2
      lambda4 <- lambda3; p[[4]] <- p[[3]]
      lambda3 <- lambda2; p[[3]] <- p[[2]]
      x1 <- log10(lambda1) 
      x4 <- log10(lambda4)
      x2 <- (x4 + phi * x1) / (1+phi) 
      x3 <- x1 + x4 - x2
      lambda2 <- 10^(x2)
      p[[2]] <- l_curve_P(glmnet_obj, y, x, weights, lambda2)
    }else{
      lambda <- lambda3
      lambda1 <- lambda2; p[[1]] <- p[[2]]
      lambda2 <- lambda3; p[[2]] <- p[[3]]
      x1 <- log10(lambda1) 
      x4 <- log10(lambda4)
      x2 <- (x4 + phi * x1) / (1+phi) 
      x3 <- x1 + x4 - x2
      lambda3 <- 10^(x3)
      p[[3]] <- l_curve_P(glmnet_obj, y, x, weights, lambda3)
    }
  }
  return(lambda)
}

##########################################################
# adaptive lasso for real-number data and complex-number data


ASG_build_library <- function(u, dx, dt, sg_para=NULL, deg_max=6, width_max=10, noise=0.01, diff_o=3, poly_o=3,
                              ref_data='Gblur', ...){
  u <- as.matrix(u)
  if(noise > 0){
    un <- u + rnorm(nrow(u)*ncol(u), 0, noise*sd(u))
  }else{
    un <- u
  }
  if(ref_data == 'Gblur'){
    un_gblur <- cv$GaussianBlur(un, as.integer(c(3,3)), 0)
  }else{
    un_gblur <- un
  }
  
  deg_min <- diff_o
  #
  n=nrow(u);m=ncol(u)
  if(is.null(sg_para)){
    deg=deg_min:deg_max#;width_max=10
    mse_t = mse_x <- matrix(Inf,nrow=deg_max,ncol=width_max)
    for(d in deg){
      width_min <- (d+3-(d %% 2)-1)/2
      for(w in width_min:width_max){
        un_sg <- matrix(0,nrow=n,ncol=m)
        for(i in 1:n){
          un_sg[i,]<-sgolayfilt(un[i,],p=d,n=2*w+1)
        }
        mse_t[d,w] <- sum((un_sg-un_gblur)^2)/(nrow(un_sg)*ncol(un_sg))
        un_sg <- matrix(0,nrow=n,ncol=m)
        for(i in 1:m){
          un_sg[,i]<-sgolayfilt(un[,i],p=d,n=2*w+1)
        }
        mse_x[d,w] <- sum((un_sg-un_gblur)^2)/(nrow(un_sg)*ncol(un_sg))
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
  
  linear_sys <- build_linear_system(un_gblur, dt, dx, D=diff_o, P=poly_o, time_diff = 'poly', space_diff = 'poly',
                                    deg_x=deg_x, deg_t=deg_t,width_x=width_x, width_t=width_t, ...)
  return(linear_sys)
}

ada_lasso_pareto_f <- function(x, y, ada_weights_max = 5, crit = 'AIC', ada_weights_gamma = NULL, ...){
  ## remove intercept when do glmnet
  # ridge_lambda = NULL; nfolders =5
  x_colnames <- colnames(x)
  n <- nrow(x)
  
  if(is.null(x_colnames)){
    x_colnames <- paste0('X',1:ncol(x))
  }
  x <- as.matrix(x)
  if(sum(x[,1]) == nrow(x)){
    x <- cbind(rnorm(n,1,1e-15), x[,-1])
  }else{
    x <- cbind(rnorm(n,1,1e-15), x)
  }
  
  if(is.null(ada_weights_gamma)){
    ada_weights_gamma_grid <- 1:ada_weights_max
  }else{
    ada_weights_gamma_grid <- ada_weights_gamma
  }
  resuslts_list <- list()
  for(g in ada_weights_gamma_grid){
    gamma_vec = bic = mse = lambda_all = numeric()
    rmse_model = active_index_list_temp = rmse_model_temp = glmnet_model = active_index_list_inter = list()
    active_index_list <- list()
    active_index_0 <- 1:ncol(x)
    ## the adaptive lasso with intercept
    active_index_list[[1]] = active_index_list_inter[[1]] <- active_index_0
    penal.coef <- lsfit(x = scale(as.matrix(x)[,-1]), y = scale(y), intercept=T, tolerance = 1e-20)[[1]] # ols fit with standardized design matrix
    weights <- 1/abs(penal.coef)^g
    model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F)
    model_lambda <- model$lambda
    # new lambda grid 0.0001 is the default, 0.0001 is the default lambda.min.ratio	
    # lambda_grid <- exp(seq(log(max(model_lambda)), log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5, length=100))
    ## try to find the lambda_max for the adaptive lasso with glmnet
    lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
    model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
    while(any(coef(model, s = model$lambda[1])[-1] != 0) ){ # glmnet sometimes has a bug, we need to approximately find the lambda_max 
      model_lambda <- model_lambda*2
      lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
      model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
    }
    glmnet_model[[1]] <- model
    len_all_model <- length(apply(model$beta, 2, function(x) length(which(x!=0))))
    if(len_all_model!=100){break}
    k1 <- 100; k4 <- 1
    while(all(coef(model,s=model$lambda[k4])[-1]==0)){
      k4=k4+1
    }
    p_lambda <- NULL
    while(is.null(p_lambda)){
      p_lambda <- tryCatch(p_lambda <- l_curve_corner_search(model,y,x,weights,model$lambda[k1],model$lambda[k4]), 
                           error = function(e) NULL)
      k1 <- k1 - 1
    }
    lambda_all[1] = p_lambda <- model$lambda[which.min(abs(p_lambda-model$lambda))]
    active_index_list_inter[[2]] <- which(coef(model,s=p_lambda)[-1]!=0)
    active_index_list[[2]] = active_index_list_inter[[2]]
    ## iteratively trimming the design matrix and implememnt the adaptive lasso with ols weights until converge
    i <- 2
    while(length(active_index_list[[i]]) != length(active_index_list[[i-1]]) & length(active_index_list[[i]])>1){
      # print(i)
      active_index <- active_index_list[[i]]
      penal.coef <- lsfit(x = scale(as.matrix(x)[,active_index]), y = scale(y), intercept=F, tolerance = 1e-20)[[1]]
      active_index_temp = active_index_list_temp = rmse_model_temp = list()
      
      # ada_weights_gamma = log10(kappa(scale(as.matrix(x)[,active_index])))
      weights <- 1/(abs(penal.coef)^g)
      model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F)
      model_lambda <- model$lambda
      lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
      model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
      while(any(coef(model, s = model$lambda[1])[-1] != 0) ){ # glmnet sometimes has a bug, we need to approximately find the lambda_max 
        model_lambda <- model_lambda*2
        lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
        model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
      }
      glmnet_model[[i]] <- model
      len_all_model <- length(apply(model$beta, 2, function(x) length(which(x!=0))))
      if(len_all_model!=100){break}
      k1 <- 100; k4 <- 1
      while(all(coef(model,s=model$lambda[k4])[-1]==0)){
        k4=k4+1
      }
      p_lambda <- NULL
      while(is.null(p_lambda)){
        p_lambda <- tryCatch(p_lambda <- l_curve_corner_search(model,y,x[,active_index],weights,model$lambda[k1],model$lambda[k4]), 
                             error = function(e) NULL)
        k1 <- k1 - 1
      }
      lambda_all[i] = p_lambda <- model$lambda[which.min(abs(p_lambda-model$lambda))]
      active_index_list_inter[[i+1]] <- active_index_list_inter[[i]][which(coef(model,s=p_lambda)[-1]!=0)]
      active_index_list[[i+1]] = active_index_list_inter[[i+1]]
      i <- i+1
    }
    active_index_list <- active_index_list[-1]
    active_index_list_length <- sapply(active_index_list, length)
    active_index_list <- active_index_list[which(active_index_list_length > 0)]
    
    crits_all <- sapply(seq_along(active_index_list), function(i){
      model_index <- which.min(abs(lambda_all[i]- glmnet_model[[i]]$lambda))
      tLL <- glmnet_model[[i]]$nulldev - deviance(glmnet_model[[i]])[model_index]
      k = df <- length(active_index_list[[i]])
      AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
      AIC <- 2*k - tLL
      BIC <- log(n)*k - tLL
      model_coef0 <- coef(glmnet_model[[i]], s=lambda_all[i])[-1]
      model_coef <- model_coef0[which(model_coef0!=0)]
      yhat <- as.matrix(x)[,active_index_list[[i]]] %*% as.matrix(model_coef)
      mse <- mean((y-yhat)^2)
      bic_zou <- mse + ((log(n) / n) * df)
      aic_zou <- mse + ((2 / n) * df)
      return(c(AIC, AICc, BIC, aic_zou, bic_zou))
    })
    crits <- switch (crit,
                     AIC = crits_all[1,],
                     AICc = crits_all[2,],
                     BIC = crits_all[3,],
                     aic_zou = crits_all[4,],
                     bic_zou = crits_all[5,],
    )
    # min_rmse_index <- which.min(sapply(rmse_model, min))+1 # the first one is the full model
    active_index_out <- active_index_list[[which.min(crits)]]
    resuslts_list[[g]] <- list(active = active_index_out, crit = min(crits))
  }
  res_crit = which.min(sapply(resuslts_list, function(e) e[[2]]))
  active_index_out <- resuslts_list[[res_crit]]$active
  
  full_coeff <- numeric(ncol(x))
  ols_coeff <- lsfit(x = as.matrix(x)[,active_index_out], y = y, intercept=F, tolerance = 1e-20)[[1]]
  full_coeff[active_index_out] <- ols_coeff
  names(full_coeff) <- x_colnames
  output <- list(full_coeff = full_coeff, active_coeff = ols_coeff, active_index_list = resuslts_list, gamma = res_crit)
  return(output)
}

ada_lasso_pareto_f_complex <- function(x, y, ada_weights_max = 5, crit = 'AIC', ada_weights_gamma = NULL, ...){
  ## remove intercept when do glmnet
  x_colnames <- colnames(x)
  n <- nrow(x)
  
  if(is.null(x_colnames)){
    x_colnames <- paste0('X',1:ncol(x))
  }
  x <- as.matrix(x)
  
  if(is.null(ada_weights_gamma)){
    ada_weights_gamma_grid <- 1:ada_weights_max
  }else{
    ada_weights_gamma_grid <- ada_weights_gamma
  }
  resuslts_list <- list()
  for(g in ada_weights_gamma_grid){
    gamma_vec = bic = mse = lambda_all = numeric()
    rmse_model = active_index_list_temp = rmse_model_temp = glmnet_model = active_index_list_inter = list()
    active_index_list <- list()
    active_index_0 <- 1:ncol(x)
    ## the adaptive lasso with intercept
    active_index_list[[1]] = active_index_list_inter[[1]] <- active_index_0
    penal.coef <- lsfit(x = scale(as.matrix(x)[,-1]), y = scale(y), intercept=T, tolerance = 1e-20)[[1]] # ols fit with standardized design matrix
    weights <- 1/abs(penal.coef)^g
    model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F)
    model_lambda <- model$lambda
    # new lambda grid 0.0001 is the default, 0.0001 is the default lambda.min.ratio	
    # lambda_grid <- exp(seq(log(max(model_lambda)), log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5, length=100))
    lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
    model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
    while(any(coef(model, s = model$lambda[1])[-1] != 0) ){
      model_lambda <- model_lambda*2
      lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
      model <- glmnet::glmnet(as.matrix(x), y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid, maxit = 1e6)
    }
    glmnet_model[[1]] <- model
    len_all_model <- length(apply(model$beta, 2, function(x) length(which(x!=0))))
    if(len_all_model!=100){break}
    # apply(model$beta, 2, function(x) length(which(x!=0)))
    k1 <- 100; k4 <- 1
    while(all(coef(model,s=model$lambda[k4])[-1]==0)){
      k4=k4+1
    }
    p_lambda <- NULL
    while(is.null(p_lambda)){
      p_lambda <- tryCatch(p_lambda <- l_curve_corner_search(model,y,x,weights,model$lambda[k1],model$lambda[k4]), 
                           error = function(e) NULL)
      k1 <- k1 - 1
      # k4 <- k4 + 1
    }
    lambda_all[1] = p_lambda <- model$lambda[which.min(abs(p_lambda-model$lambda))]
    active_index_list_inter[[2]] <- which(coef(model,s=p_lambda)[-1]!=0)
    active_index_list[[2]] = active_index_list_inter[[2]]
    i <- 2
    while(length(active_index_list[[i]]) != length(active_index_list[[i-1]]) & length(active_index_list[[i]])>1){
      # print(i)
      active_index <- active_index_list[[i]]
      penal.coef <- lsfit(x = scale(as.matrix(x)[,active_index]), y = scale(y), intercept=F, tolerance = 1e-20)[[1]]
      active_index_temp = active_index_list_temp = rmse_model_temp = list()
      
      weights <- 1/(abs(penal.coef)^g)
      model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F)
      model_lambda <- model$lambda
      lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
      model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid)
      while(any(coef(model, s = model$lambda[1])[-1] != 0) ){
        model_lambda <- model_lambda*2
        lambda_grid <- exp(seq(log(max(model_lambda)), max(log(1e-7),log(max(model_lambda)*0.0001)-log(max(model_lambda))*1.5), length=100))
        model <- glmnet::glmnet(as.matrix(x)[,active_index], y, penalty.factor = weights, alpha=1,intercept = F, lambda = lambda_grid, maxit = 1e6)
      }
      glmnet_model[[i]] <- model
      len_all_model <- length(apply(model$beta, 2, function(x) length(which(x!=0))))
      if(len_all_model!=100){break}
      k1 <- 100; k4 <- 1
      while(all(coef(model,s=model$lambda[k4])[-1]==0)){
        k4=k4+1
      }
      p_lambda <- NULL
      while(is.null(p_lambda)){
        p_lambda <- tryCatch(p_lambda <- l_curve_corner_search(model,y,x[,active_index],weights,model$lambda[k1],model$lambda[k4]), 
                             error = function(e) NULL)
        k1 <- k1 - 1
      }
      lambda_all[i] = p_lambda <- model$lambda[which.min(abs(p_lambda-model$lambda))]
      active_index_list_inter[[i+1]] <- active_index_list_inter[[i]][which(coef(model,s=p_lambda)[-1]!=0)]
      active_index_list[[i+1]] = active_index_list_inter[[i+1]]
      i <- i+1
    }
    active_index_list <- active_index_list[-1]
    active_index_list_length <- sapply(active_index_list, length)
    active_index_list <- active_index_list[which(active_index_list_length > 0)]
    
    crits_all <- sapply(seq_along(active_index_list), function(i){
      model_index <- which.min(abs(lambda_all[i]- glmnet_model[[i]]$lambda))
      tLL <- glmnet_model[[i]]$nulldev - deviance(glmnet_model[[i]])[model_index]
      k = df <- length(active_index_list[[i]])
      AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
      AIC <- 2*k - tLL
      BIC <- log(n)*k - tLL
      model_coef0 <- coef(glmnet_model[[i]], s=lambda_all[i])[-1]
      model_coef <- model_coef0[which(model_coef0!=0)]
      yhat <- as.matrix(x)[,active_index_list[[i]]] %*% as.matrix(model_coef)
      mse <- mean((y-yhat)^2)
      bic_zou <- mse + ((log(n) / n) * df)
      aic_zou <- mse + ((2 / n) * df)
      return(c(AIC, AICc, BIC, aic_zou, bic_zou))
    })
    crits <- switch (crit,
                     AIC = crits_all[1,],
                     AICc = crits_all[2,],
                     BIC = crits_all[3,],
                     aic_zou = crits_all[4,],
                     bic_zou = crits_all[5,],
    )
    # min_rmse_index <- which.min(sapply(rmse_model, min))+1 # the first one is the full model
    active_index_out <- active_index_list[[which.min(crits)]]
    resuslts_list[[g]] <- list(active = active_index_out, crit = min(crits))
  }
  res_crit = which.min(sapply(resuslts_list, function(e) e[[2]]))
  active_index_out <- resuslts_list[[res_crit]]$active
  
  full_coeff <- numeric(ncol(x))
  ols_coeff <- lsfit(x = as.matrix(x)[,active_index_out], y = y, intercept=F, tolerance = 1e-20)[[1]]
  full_coeff[active_index_out] <- ols_coeff
  names(full_coeff) <- x_colnames
  output <- list(full_coeff = full_coeff, active_coeff = ols_coeff, active_index_list = resuslts_list, gamma = res_crit)
  return(output)
}
