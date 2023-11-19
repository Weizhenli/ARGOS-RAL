## needed packages ----------------------
library(pracma)
library(plotly)
library(R.matlab)
library(glmnet)
library(dplyr)
library(leaps)
library(signal)
library(orthopolynom)
library(reticulate)
library(signal)
# STRidge ---------------------------
STRidge <- function(X0, y, lam, maxit, tol, normalize = 2, print_results = F) {
  "
  Sequential Threshold Ridge Regression algorithm for finding (hopefully) sparse 
  approximation to X^{-1}y.  The idea is that this may do better with correlated observables.
  This assumes y is only one column
  "
  # for debug
  #X0=R; y=Ut; lam=lam; maxit=25; tol=11; normalize = 2; print_results = F
  n <- dim(X0)[1]
  d <- dim(X0)[2]
  if(!is.complex(as.matrix(X0)[1,1])){
    X <- matrix(0, nrow = n, ncol = d)
    num_type='real'
  }else{
    X <- matrix(complex(), nrow = n, ncol = d)
    y <- as.complex(y)
    num_type='complex'
  }
  #
  # First normalize data
  if(normalize != 0){
    Mreg = numeric(d)
    for(i in 1:d){
      Mreg[i] <- 1/(pracma::Norm(X0[,i],normalize))
      X[,i] <- Mreg[i]*X0[,i]
    }
  }else{X=X0}
  
  # Get the standard ridge estimate
  if(lam != 0){
    y.temp1 <- t(X)%*%y
    X.temp1 <- t(X)%*%X + lam*diag(d)
    #w<-solve(X.temp1,y.temp1)
    #w <- solve(t(X.temp1)%*%X.temp1)%*%t(X.temp1)%*%y.temp1
    #w <- solve(Conj(t(X.temp1))%*%X.temp1)%*%Conj(t(X.temp1))%*%y.temp1
    w <- lsfit(X.temp1,y.temp1,intercept=F,tolerance = 1e-20)[[1]]
    #w = lm(y.temp1~ 0 + X.temp1)$coef
  }else{
    w <- lsfit(X,y,intercept=F,tolerance = 1e-20)[[1]]
    #w <- solve(Conj(t(X)) %*% X) %*% Conj(t(X)) %*% y
    #w <- lm(y~ 0 + X)$coef
  }
  num_relevant <- d
  biginds <- which(abs(w) > tol)
  
  # Threshold and continue
  for(j in 1:maxit){
    #j=1
    # Figure out which items to cut out
    smallinds <- which(abs(w) < tol)
    new_biginds <- (1:d)[-smallinds]
    
    # If nothing changes then stop
    if(num_relevant == length(new_biginds)){
      break
    }else{num_relevant <- length(new_biginds)}
    #print(j)
    # Also make sure we didn't just lose all the coefficients
    if(length(new_biginds) == 0){
      if(j == 1){
        #if print_results: print "Tolerance too high - all coefficients set below tolerance"
        #w
        return(w)
      }else{break}
    }
    biginds <- new_biginds
    #print(j)
    # Otherwise get a new guess
    w[smallinds] <- 0
    if(lam != 0){
      y.temp2 <- t(X[,biginds])%*%y
      X.temp2 <- t(X[,biginds])%*%X[,biginds] + lam*diag(length(biginds))
      w[biginds] <- lsfit(X.temp2,y.temp2,intercept=F,tolerance = 1e-20)[[1]]
      #w[biginds] <- solve(Conj(t(X.temp2)) %*% X.temp2) %*% Conj(t(X.temp2)) %*% y.temp2
      #w[biginds] <- lm(y.temp2~ 0 + X.temp2)$coef
    }else{
      w[biginds] <- lsfit(X[,biginds],y,intercept=F,tolerance = 1e-20)[[1]]
      #w[biginds] <- solve(Conj(t(X[,biginds])) %*% X[,biginds]) %*% Conj(t(X[,biginds])) %*% y
      #w[biginds] <- lm(y~ 0 + X[,biginds])$coef
    }
  }
  
  # Now that we have the sparsity pattern, use standard least squares to get w
  if(length(biginds)!=0){
    w[biginds] <-lsfit(X[,biginds],y,intercept=F,tolerance = 1e-20)[[1]]
    #w[biginds] <- solve(Conj(t(X[,biginds])) %*% X[,biginds]) %*% Conj(t(X[,biginds])) %*% y
    #w[biginds] <- lm(y~ 0 + X[,biginds])$coef
  }
  if(normalize != 0){
    #Mreg * w
    return(Mreg * w)
  }else{
    #w
    return(w)
  }
}
#
# TrainSTRidge ------------------------------------
TrainSTRidge <- function(R, Ut, lam, d_tol, maxit = 25, STR_iters = 10, l0_penalty = NA, 
                         normalize = 2, split = 0.8, print_best_tol = F){
  "
  This function trains a predictor using STRidge.
  
  It runs over different values of tolerance and trains predictors on a training set, then evaluates them 
  using a loss function on a holdout set.
  
  Please note published article has typo.  Loss function used here for model selection evaluates fidelity using 2-norm,
  not squared 2-norm.
  "
  # Split data into 80% training and 20% test, then search for the best tolderance.
  
  # for debug
  #R=R; Ut=Ut; maxit = 25; STR_iters = 10; l0_penalty = NA; normalize = 2; split = 0.8; print_best_tol = F
  
  # set.seed(1) # for consistency
  R.name <- names(R)
  R <- as.matrix(R)
  Ut <- as.matrix(Ut)
  n <- dim(R)[1]
  n.ind <- 1:n
  train <- sample(n.ind, n*split) # index of training data random sampling
  test <- n.ind[-train]
  TrainR <- R[train,]
  TestR <- R[test,]
  TrainY <- Ut[train,]
  TestY <- Ut[test,]
  D <- dim(TrainR)[2]  
  
  # Set up the initial tolerance and l0 penalty
  tol <- d_tol
  if(is.na(l0_penalty)){
    #l0_penalty <- 0.001*pracma::cond(as.matrix(R))
    # s <- svd(R)$d
    l0_penalty <- 0.001*kappa(R)
    #l0_penalty <- 0.001*np$linalg$cond(R)
  }
  
  # Get the standard least squares estimator
  w <- numeric(D)
  note <- t(TrainR)%*%TrainR
  invers.note <- tryCatch(solve(note), error=function(e) NULL)
  #note <- solve(t(TrainR)%*%TrainR)%*%t(TrainR)%*%TrainY
  #note <- Matrix::solve(TrainR,TrainY)
  if(is.null(invers.note)){
    #w_best <- np$linalg$lstsq(TrainR, TrainY)[[1]]
    TrainR <- TrainR+0i
    TrainY <- TrainY+0i
    # w_best <- solve(Conj(t(TrainR))%*%TrainR)%*%Conj(t(TrainR))%*%TrainY
    # w_best <- Re(w_best)
    w_best <- lsfit(TrainR, TrainY, intercept=F, tolerance = 1e-16)[[1]]
  }else{w_best <- solve(note)%*%t(TrainR)%*%TrainY}
  
  #w_best <- solve(Conj(t(TrainR))%*%TrainR)%*%Conj(t(TrainR))%*%TrainY
  #w_best <- lm(TrainY~ 0 + TrainR)$coef
  #w_best <- np$linalg$lstsq(TrainR, TrainY)[[1]] # use python least square python deals with singular fit
  n_none_0NA_w_best <- length(which(Mod(w_best)> 1e-13 | Mod(w_best)< (-1e-13))) # 0 in lm()$coef is not exact 0. 
  #err_best <- pracma::Norm(TestY - TestR%*%w_best)^2 + l0_penalty * n_none_0NA_w_best
  err_best <- pracma::Norm(TestY - TestR%*%w_best) + l0_penalty * n_none_0NA_w_best
  tol_best = 0
  
  # Now increase tolerance until test performance decreases
  for(iter in 1:maxit){
    #iter=10
    # Get a set of coefficients and error
    #print(iter)
    w <- STRidge(R,Ut,lam,STR_iters,tol,normalize = normalize)
    n_none_0NA_w <- length(which(Mod(w)>1e-13 | Mod(w)<(-1e-13)))
    #err <- pracma::Norm(TestY - TestR%*%matrix(w,ncol=1))^2 + l0_penalty*n_none_0NA_w # can use other way for errors
    err <- pracma::Norm(TestY - TestR%*%matrix(w,ncol=1)) + l0_penalty*n_none_0NA_w # can use other way for errors
    
    # Has the accuracy improved?
    if(err <= err_best){
      err_best <- err
      w_best <- w
      tol_best <- tol
      tol <- tol + d_tol
    }else{
      tol <- max(0,tol-2*d_tol)
      d_tol <- 2*d_tol / (maxit - iter)
      tol <- tol + d_tol
    }
  }
  if(print_best_tol){
    print(paste("Optimal tolerance:", tol_best))
  }
  return(w_best)
}
TrainSTRidge_boot <- function(data,lam, d_tol, ...){
  TrainSTRidge(R=data[,-1], Ut=data[,1], lam, d_tol, ...)
}


# ConvSmoother --------------------------
ConvSmoother <- function(x, p, sigma){
  "
  Smoother for noisy data
  Inpute = x, p, sigma
  x = one dimensional series to be smoothed
  p = width of smoother
    sigma = standard deviation of gaussian smoothing kernel
  "
  n <- length(x)
  y <- numeric(n)
  g <- exp(-seq(-p, p, length.out = 2*p)^2/(2.0*sigma^2))
  
  for(i in 1:n){
    a <- max(i-p, 1) # python from 0, r from 1
    b <- min(i-1+p, n)
    c <- max(1, p-i+2)
    d <- min(2*p, p+n-i+1)
    #print(c(a,b,c,d))
    y[i] <- sum(x[a:b] * g[c:d])/sum(g[c:d])
  }
  return(y)
}

#ConvSmoother(1:10,6,2)

# TikhonovDiff -----------------------------------
TikhonovDiff <- function(f, dx, lam, d = 1){
  
  #Tikhonov differentiation.
  #return argmin_g \|Ag-f\|_2^2 + lam*\|Dg\|_2^2
  #where A is trapezoidal integration and D is finite differences for first dervative
  #It looks like it will work well and does for the ODE case but 
  #tends to introduce too much bias to work well for PDEs.  If the data is noisy, try using
  #polynomials instead.
  
  # Initialize a few things
  n <- length(f)
  if(!is.complex(f[1])){
    num_type='real'
  }else{num_type='complex'}
  f <- f - f[1]
  
  # Get a trapezoidal approximation to an integral
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 2:n){
    A[i,i] <- dx/2
    A[i,0] <- dx/2
    for(j in 2:i){
      A[i,j] <- dx
    }
  }
  
  #e <- rep(1,n-1)
  sparse.m <- function(n, e){
    # Further, try to construct function like scipy.sparse.diags() in Python 
    i <- rep(1:(n-1), each = 2)
    j <- c(1, rep(2:(n-1), each = 2), n)
    x <- rep(c(-e, e), times = (n-1))
    return(as.matrix(Matrix::sparseMatrix(i, j, x = x)))
  }
  D <- sparse.m(n, 1)
  
  # Invert to find derivative
  y.temp <- t(A)%*%f #A.T.dot(f)
  x.temp <- t(A)%*%A + lam*t(D)%*%D #A.T.dot(A) + lam*D.T.dot(D)
  g <- lsfit(x.temp,y.temp,intercept=F,tolerance = 1e-20)[[1]]
  #g <-solve(t(x.temp) %*% x.temp) %*% t(x.temp) %*% y.temp
  #g <- lm(y.temp~ 0 + x.temp)$coef
  if(d == 1){
    #g
    return(g)
  }else{
    # If looking for a higher order derivative, this one should be smooth so now we can use finite differences
    return(FiniteDiff(g, dx, d-1))
  }
}
#TikhonovDiff(1:10,6,2)

# FiniteDiff ------------------------------
FiniteDiff <- function(u, dx, d){
  "
  Takes dth derivative data using 2nd order finite difference method (up to d=3)
  Works but with poor accuracy for d > 3
  Input:
  u = data to be differentiated
  dx = Grid spacing.  Assumes uniform spacing
  "
  # for debug
  #u <-u[1,];dx;d=1
  
  u<-as.matrix(u)
  if(!is.complex(u[1,1])){
    u <- as.numeric(u)
  }else{
    u <- as.complex(u) # important
  }
  #
  FiniteDiff_123 <- function(u, dx, d){
    n <- length(u)
    u<-as.matrix(u)
    if(!is.complex(u[1,1])){
      ux <- numeric(n)
    }else{
      ux <- complex(n)
    }
    #
    if(d == 1){
      for(i in 2:(n-1)){
        ux[i] <- (u[i+1]-u[i-1]) / (2*dx)
      }
      ux[1] <- (-3.0/2*u[1] + 2*u[2] - u[3]/2) / dx
      ux[n] <- (3.0/2*u[n] - 2*u[n-1] + u[n-2]/2) / dx
      return(ux)
    }
    if(d == 2){
      for(i in 2:(n-1)){
        ux[i] <- (u[i+1]-2*u[i]+u[i-1]) / dx^2
      }
      ux[1] <- (2*u[1] - 5*u[2] + 4*u[3] - u[4]) / dx^2
      ux[n] <- (2*u[n] - 5*u[n-1] + 4*u[n-2] - u[n-3]) / dx^2
      return(ux)
    }
    if(d == 3){
      for(i in 3:(n-2)){
        ux[i] = (u[i+2]/2-u[i+1]+u[i-1]-u[i-2]/2) / dx^3
      }
      ux[1] = (-2.5*u[1]+9*u[2]-12*u[3]+7*u[4]-1.5*u[5]) / dx^3
      ux[2] = (-2.5*u[2]+9*u[3]-12*u[4]+7*u[5]-1.5*u[6]) / dx^3
      ux[n] = (2.5*u[n]-9*u[n-1]+12*u[n-2]-7*u[n-3]+1.5*u[n-4]) / dx^3
      ux[n-1] = (2.5*u[n-1]-9*u[n-2]+12*u[n-3]-7*u[n-4]+1.5*u[n-5]) / dx^3
      return(ux)
    }
  }
  temp<-c()
  y <- u
  temp1 <- u
  while(d > 3){
    temp <- FiniteDiff_123(temp1,dx,3)
    d <- d-3
    temp1 <- temp
  }
  while(d == 3){
    if(!is.null(temp)){
      y <- temp
    }
    return(FiniteDiff_123(y,dx,3))
    break
  }
  while(d == 2){
    if(!is.null(temp)){
      y <- temp
    }
    return(FiniteDiff_123(y,dx,2))
    break
  }
  while(d == 1){
    if(!is.null(temp)){
      y <- temp
    }
    return(FiniteDiff_123(y,dx,1))
    break
  }
}

#set.seed(1)
#u=rnorm(10)
# test
#for(i in 1:9){
#  print(FiniteDiff(u,1,i))
#}

# PolyDiff ----------------------
# a python function is needed
PolyDiff <- function(y, dx, deg=3, width=5, diff=1){
  #dx=diff(x)[1];deg=5;width=5;diff=1
  n <- length(y)
  if(diff!=0){
    sgdy <- matrix(NA,nrow=n-2*width,ncol=diff)
    for(d in 1:diff){
      #sgolayfilt
      sgdy[,d] <- ((factorial(d)*pracma::savgol(y, width*2+1, deg, d))/dx^d)[(width+1):(n-width)]
    }
  }else{
    sgdy <- (pracma::savgol(y, width*2+1, deg))[(width+1):(n-width)]
  }
  return(sgdy)
}
#=seq(1,10,length.out = 20)
#u=seq(1,10,length.out = 30)*4
#drt=c('D:/Users/Administrator/anaconda3/envs/python37')
#PolyDiff(u,x,drt=drt)

# PolyDiffPoint -------------------
PolyDiffPoint <- function(u, x, deg = 3, dif = 1, index = NA){
  #library(reticulate)
  #use_python(drt)
  #np<-import('numpy')
  "
  Same as above but now just looking at a single point

  u = values of some function
  x = x-coordinates where values are known
  deg = degree of polynomial to use
  diff = maximum order derivative we want
  "
  #u=W[x+1,y+1,(t-(N-1)/2+1):(t+(N+1)/2)];x=(1:N-1)*dt;dif=1;index = NA
  n <- length(x)
  degint <- as.integer(deg)
  if(is.na(index)){index <- as.integer((n-1)/2)}
  poly <- np$polynomial$chebyshev$Chebyshev$fit(x,u,degint)
  # Fit to a Chebyshev polynomial
  # better conditioned than normal polynomials
  
  # Take derivatives
  deriva <- numeric()
  for(d in 2:(dif+1)){
    dint <- as.integer(d-1)
    deriva[d] <- np$polynomial$chebyshev$Chebyshev$deriv(poly,m=dint)(x[index+1])
  }
  return(deriva[-1])
}

library(orthopolynom)
PolyDiffPoint <- function(y, x, deg = 3, dif = 1, index = NA, type='c'){
  n <- length(x)
  if(is.na(index)){index <- as.integer((n-1)/2)}
  
  if(type=='c'){
    cheb <- orthopolynom::chebyshev.c.polynomials(deg)
  }else if(type=='s'){
    cheb <- orthopolynom::chebyshev.s.polynomials(deg)
  }else if(type=='t'){
    cheb <- orthopolynom::chebyshev.t.polynomials(deg)
  }else if(type=='u'){
    cheb <- orthopolynom::chebyshev.u.polynomials(deg)
  }
  
  cheb.x <- cbind(rep(1,length(x)),sapply(2:(deg+1),function(ls) {eval(str2expression(paste(cheb[[ls]])))})) # evaluate cheb.poly with x
  poly <- lsfit(cheb.x, y, intercept=F,tolerance=1e-20)[[1]] # find the poly coeffs and solve singular matrix
  
  # Take derivatives
  deriva <- numeric()
  for(d in 2:(dif+1)){
    deriva[d] <- poly_deriv(cheb,poly,d-1,x[index+1])
  }
  return(deriva[-1])
}


D.fun <- function(cheb,item,k=1){ # find the kth derivative for the chebyshev polynomial object 
  D.expr <- D(str2expression(paste(cheb[[item]])),'x')
  i <- 1
  while(k>i){
    D.expr <- D(D.expr,'x')
    i <- i+1
  }
  return(D.expr)
}

PolyDiffPoint2 <- function(y, x, deg = 3, dif = 1, index = NA){
  #deg=4;dif=2
  n <- length(x)
  if(is.na(index)){index <- as.integer((n-1)/2)}
  cheb <- chebyshev.c.polynomials(deg)
  cheb.x <- cbind(rep(1,length(x)),sapply(2:(deg+1),function(ls) eval(str2expression(paste(cheb[[ls]])))))
  poly <- solve(t(cheb.x)%*%cheb.x)%*%t(cheb.x)%*%y
  #cheb.x[index+1,]%*%poly
  D.fun <- function(item,k=1){ # find the kth derivative of 
    D.expr <- D(str2expression(paste(cheb[[item]])),'x')
    i <- 1
    while(k>i){
      D.expr <- D(D.expr,'x')
      i <- i+1
    }
    return(D.expr)
  }
  #D.fun(5,2)
  #D.beta <- lapply(1:(deg+1),function(ls) D.fun(ls,dif))
  #D.beta <- lapply(1:(deg+1),function(ls) D(str2expression(paste(cheb[[ls]])),'x'))
  eval_D.poly <- function(x,k=1){
    D.poly <- lapply(1:(deg+1),function(ls) D.fun(ls,k))
    sum(sapply(1:(deg+1), function(ls) eval(D.poly[[ls]])*poly[ls]))
    #(eval(D.beta[[1]])*poly[1]+eval(D.beta[[2]])*poly[2]+eval(D.beta[[3]])*poly[3]+eval(D.beta[[4]])*poly[4]+eval(D.beta[[5]])*poly[5])[index+1]
  }
  deriva <- numeric()
  for(d in 1:dif){
    deriva[d] <- eval_D.poly(x[index+1],k=d)
  }
  return(deriva)
}
#x=seq(1,10,length.out = 20)
#u=seq(1,10,length.out = 20)*4
#PolyDiffPoint(u,x,index=1,drt=drt)

## Discrete Fourier Transform sample frequencies --------------
fftfreq <- function(n,d=1){
  # calculate the Discrete Fourier Transform sample frequencies
  val = 1.0 / (n * d)
  results = numeric(n)
  N = (n-1) %/% 2 + 1
  p1 = seq(0, N-1)
  results[1:N] = p1
  p2 = seq(-(n %/% 2), 0-1)
  results[(N+1):length(c(p1,p2))] = p2
  results * val
}

#np$fft$fftfreq(as.integer(7))==fftfreq(7)

# build_linear_system --------------------
build_linear_system <- function(u, dt, dx, D = 3, P = 3, Dt=1, time_diff = 'poly', space_diff = 'poly',
                                lam_t = NA, lam_x = NA, width_x = NA, width_t = NA, 
                                deg_x = 5, deg_t = NA, sigma = 2){
  u <- as.matrix(u) # u must be matrix
  if(!is.complex(u[2,1])){
    number_type <- 'real'
    #ut <- matrix(0, nrow = n2, ncol = m2)
  }else{
    number_type <- 'complex'
    #ut <- matrix(complex(), nrow = n2, ncol = m2)
  }
  
  
  n <- dim(u)[1]
  m <- dim(u)[2]
  if(is.na(width_x)){width_x <- n/10}
  if(is.na(width_t)){width_t <- m/10}
  if(is.na(deg_t)){deg_t <- deg_x}
  
  # If we're using polynomials to take derviatives, then we toss the data around the edges.
  if(time_diff == 'poly'){
    m2 <- m-2*width_t
    offset_t <- width_t + 1
  }else{
    m2 <- m
    offset_t <- 1
  }
  if(space_diff == 'poly'){
    n2 <- n-2*width_x
    offset_x <- width_x + 1
  }else{
    n2 <- n
    offset_x <- 1
  }
  if(is.na(lam_t)){lam_t <- 1/m}
  if(is.na(lam_x)){lam_x <- 1/n}
  
  ########################
  # First take the time derivaitve for the left hand side of the equation
  ########################
  
  if(number_type=='real'){
    ut <- matrix(0, nrow = n2, ncol = m2)
  }else{
    ut <-  matrix(complex(), nrow = n2, ncol = m2)
  }
  #
  if(time_diff == 'FDconv'){
    if(number_type=='real'){
      Usmooth <- matrix(0, nrow = n, ncol = m)
    }else{
      Usmooth <- matrix(complex(), nrow = n, ncol = m)
    }
    #Usmooth <- matrix(0, nrow = n, ncol = m)
    # Smooth across x cross-sections
    for(j in 1:m){
      Usmooth[,j] <- as.numeric(ConvSmoother(u[,j],width_t,sigma))
    }
    # Now take finite differences
    for(i in 1:n2){
      ut[i,] <- as.numeric(FiniteDiff(Usmooth[i + offset_x - 1,],dt,Dt))
    }
  }else if(time_diff == 'poly'){
    polyx <- seq(0, (m-1)*dt, length.out=m)
    for(i in 1:n2){
      #ut[i,] <- as.numeric(PolyDiff(u[i + offset_x - 1,],polyx,dif=Dt,width=width_t,deg=deg_t)[,1])
      ut[i,] <- as.numeric(PolyDiff(u[i + offset_x - 1,],dt,deg=deg_t,width=width_t,diff=Dt)[,Dt])
    }
  }else if(time_diff == 'Tik'){
    for(i in 1:n2){
      ut[i,] = as.numeric(TikhonovDiff(u[i + offset_x - 1,], dt, lam_t, d=Dt))
    }
  }else{
    for(i in 1:n2){
      ut[i,] = as.numeric(FiniteDiff(u[i + offset_x - 1 ,],dt,Dt))
    }
  }
  #ut1
  dim(ut) <- c(n2*m2,1)

  ########################
  # Now form the rhs one column at a time, and record what each one is
  ########################
  u2 <- u[(offset_x:(n-offset_x+1)),(offset_t:(m-offset_t+1))]
  if(number_type=='real'){
    Theta <- matrix(0, nrow = n2*m2, ncol = ((D+1)*(P+1)))
    ux <- matrix(0, nrow = n2, ncol = m2)
  }else{
    Theta <- matrix(complex(), nrow = n2*m2, ncol = ((D+1)*(P+1)))
    ux <- matrix(complex(), nrow = n2, ncol = m2)
  }
  rhs_description <- character((D+1)*(P+1))
  
  if(space_diff == 'poly'){
    poly_x <- seq(0,(n-1)*dx,length.out=n)
    lenDU <- n - 2*width_x
    Du <- array(0,dim = c(lenDU,D,m2))
    for(i in 1:m2){
      Du[ , ,i] <- PolyDiff(u[,i+offset_t-1], dx, diff=D,width=width_x,deg=deg_x)
    }
  }
  if(space_diff == 'Fourier'){
    ik = 1i*fftfreq(n)*n # change to complex number
  }
  
  for(d in 1:(D+1)){
    if(d > 1){
      for(i in 1:m2){
        if(space_diff == 'Tik'){
          ux[,i] <- as.numeric(TikhonovDiff(u[,i+offset_t-1], dx, lam_x, d=(d-1)))
        }else if(space_diff == 'FDconv'){
          Usmooth <- as.numeric(ConvSmoother(u[,i+offset_t-1],width_x,sigma))
          ux[,i] <- as.numeric(FiniteDiff(Usmooth,dx,d))
        }else if(space_diff == 'FD'){
          ux[,i] <- as.numeric(FiniteDiff(u[,i+offset_t-1],dx,(d-1)))
        }else if(space_diff == 'poly'){
          ux[,i] <- as.numeric(Du[,d-1,i])
        }else if(space_diff == 'Fourier'){
          uxff <- ik^d*fft(u[,i])
          ux[,i] <- fft(uxff, inverse = T) / length(uxff)   # problem in the length of ik^d and ux[,i] change length of ik to n2
        }
      }
    }else{
      ux <- matrix(1, nrow = n2, ncol = m2)
    }
    #print(d)
    d1 <- d - 1
    for(r in 1:(P+1)){
      tempp <- ux*(u2)^(r-1)
      dim(tempp) <- c(n2*m2,1) # tempp must be in matrix frame
      Theta[, d1*(P+1)+r] <- tempp
      #Theta[:, d*(P+1)+p] = np.reshape(np.multiply(ux, np.power(u2,p)), (n2*m2), order='F')
      if(r == 2){
        rhs_description[d1*(P+1)+r] <- paste(rhs_description[d1*(P+1)+r],'u', sep = '')
      }else if(r > 2){
        rhs_description[d1*(P+1)+r] <- paste(rhs_description[d1*(P+1)+r],'u^',as.character(r-1),sep='')
      }
      if(d1 > 0){
        rhs_description[d1*(P+1)+r] = paste(rhs_description[d1*(P+1)+r], 'u_{', '', paste(rep('x', d1), collapse = ''), '}', sep = '')
      }
    }
  }
  return(list(ut=ut,Theta=Theta,rhs=rhs_description))
}



## build_Theta ---------------------------
build_Theta <- function(data, derivatives=NULL, derivatives_description=NULL, P, data_description = NULL){
  data <- as.data.frame(data)
  data.dim <- dim(data)
  n <- data.dim[1]
  d <- data.dim[2]
  if(!is.null(data_description)){
    if(length(data_description) != d){
      stop('data descrption error')
    }else{
      names(data) <- data_description
    }
  }
  if(P > 1){
    data.ploy <- as.data.frame(stats::poly(as.matrix(data),P,raw=T))
    # names of poly
    ## for the only poly
    name.poly <- character()
    for(i in names(data.ploy)){
      power <- as.character(unlist(strsplit(i, split = "")))
      name.temp <- numeric()
      for(j in 1:length(data_description)){
        if(power[2*j-1] == '0'){
          name.temp[j] <- ''
        }else if(power[2*j-1] == '1'){
          name.temp[j] <- paste(data_description[j])
        }else{name.temp[j] <- paste(data_description[j],'^',as.character(power[2*j-1]),sep='')}
      }
      name.poly[i] <- paste(name.temp,sep='',collapse = '')
    }
  }else{
    data.ploy <- as.data.frame(data)
    name.poly <- names(data.ploy)
  }

  # if has derivatives
  if(!is.null(derivatives)){
    derivatives <- as.data.frame(derivatives)
    names(derivatives) <- derivatives_description
    derivatives.dim <- dim(derivatives)
    m <- derivatives.dim[1]
    d2 <- derivatives.dim[2]
    if(n != m){stop('dimension error')}
    
    Theta <- data.frame(derivatives)
    for(i in 1:d2){
      temp <- data.ploy * derivatives[,i]
      Theta <- cbind.data.frame(Theta,temp)
    }
    
    ## add derivatives
    names.n <- numeric()
    for(i in 1:d2){
      names.j <- numeric()
      for(j in 1:length(name.poly)){
        names.j[j] <- paste(name.poly[j],derivatives_description[i],sep='')
      }
      names.n <- c(names.n,names.j)
    }
    name.poly <- c(derivatives_description,names.n)
    names(Theta) <- name.poly
    return(Theta) 
  }else{
    names(data.ploy) <- name.poly
    return(data.ploy) 
  }
}  


## reshape matrix dimension ----------------------
reshape.dim <- function(x,nrow,ncol){ # x must be a matrix or dataframe
  x.copy <- as.matrix(x)
  dim(x.copy) <- c(nrow,ncol)
  return(x.copy)
}
## only transfer data ----------------------
break1 = function(X) {
  #do.call(c, lapply(X, function(x) { c(Re(x), Im(x)) }))
  as.numeric(sapply(X, function(x) { c(Re(x), Im(x)) }))
}
break2 = function(X) {
  #do.call(c, lapply(X, function(x) { c(-Im(x), Re(x)) }))
  as.numeric(sapply(X, function(x) { c(-Im(x), Re(x)) }))
}
transfer_old <- function(Y,X){
  # Split into real variables
  
  YF = break1(Y)
  XF.List = do.call(c, lapply(apply(X, 2, as.list),
                              function(x) { list(break1(x), break2(x)) } ))
  #XF.List <- data.frame(apply(data.frame(apply(X1, 2, list)),2,function(x){list(break1(x), break2(x))}))
  #XF.List <- sapply(apply(X1, 2, as.list),function(x){list(break1(x), break2(x))})
  #system.time(XF.List <- do.call(c, lapply(apply(X1, 2, as.list),function(x) { list(break1(x), break2(x)) } )))
  #system.time(XF <- apply(X1,2,function(x){list(break1(x), break2(x))}))
  #system.time(XF2 <- sapply(apply(X1, 2, as.list),function(x){list(break1(x), break2(x))}))
  #system.time(XF2<-data.frame(apply(data.frame(apply(X1, 2, list)),2,function(x){list(break1(x), break2(x))})))
  # ...
  # Make the data.frame
  Data = data.frame(Y = YF)
  X.Names = paste('X', 1:length(XF.List), sep='')
  for (N in seq_along(XF.List)) {
    Data[[ X.Names[[N]] ]] = XF.List[[N]]
    #Data[[ X.Names[[N]] ]] = XF.List[,N]
  }
  d <- dim(X)[2]*2
  X.names <- names(X)
  coeff.names <- paste('X', 1:d, sep='')
  for(a in 1:(d/2)){
    coeff.names[2*a-1] <-  paste('Re(',X.names[a],')',sep='')
    coeff.names[2*a] <- paste('Im(',X.names[a],')',sep='')
  }
  names(Data) <- c('Y',coeff.names)
  return(Data)
}
transfer <- function(Y,X){
  # Split into real variables
  
  YF = break1(Y)
  XF.List = do.call(c, lapply(apply(X, 2, as.list),
                              function(x) { list(break1(x), break2(x)) } ))
  # Make the data.frame
  Data = data.frame(Y = YF)
  X.Names = paste('X', 1:length(XF.List), sep='')
  for (N in seq_along(XF.List)) {
    Data[[ X.Names[[N]] ]] = XF.List[[N]]
    #Data[[ X.Names[[N]] ]] = XF.List[,N]
  }
  d <- dim(X)[2]*2
  X.names <- names(X)
  coeff.names <- paste('X', 1:d, sep='')
  for(a in 1:(d/2)){
    coeff.names[2*a-1] <-  paste(X.names[a],sep='')
    coeff.names[2*a] <- paste('i',X.names[a],sep='')
  }
  names(Data) <- c('Y',coeff.names)
  return(Data)
}
## complex regression for lm() ---------------------------
fit.complex = function(Y, X) { # it can work but the coeff names and the input should be improved. 
  X.List <- apply(X, 2, as.list)
  # Split into real variables
  YF = break1(Y)
  XF.List = do.call(c, lapply(X.List,
                              function(x) { list(break1(x), break2(x)) } )) # save memory size
  # Make the data.frame
  Data = data.frame(Y = YF)
  X.Names = paste('X', 1:length(XF.List), sep='')
  X.Names2 <- paste('X', 0:(length(XF.List)/2-1), sep='')
  for (N in seq_along(XF.List)) {
    Data[[ X.Names[[N]] ]] = XF.List[[N]]
  }
  
  # Formula + Model
  Formula = paste("Y ~ ", paste(X.Names, collapse='+'), "-1")
  Model = lm(as.formula(Formula), data=Data)
  # Make them complex again
  Coeffs = sapply(seq_along(X.List),
                  function(N) {
                    ( Model$coefficients[[ X.Names[[2*N-1]] ]]
                      + Model$coefficients[[ X.Names[[2*N]] ]]*1i )
                  })
  Coeff.names <- X.Names
  if(is.null(names(X.List))){
    for(j in seq_along(X.List)){
      Coeff.names[2*j-1] <-  paste('Re(',X.Names2[j],')',sep='')
      Coeff.names[2*j] <- paste('Im(',X.Names2[j],')',sep='')
    }
  }else{
    for(j in seq_along(X.List)){
      Coeff.names[2*j-1] <-  paste('Re(',names(X.List)[j],')',sep='')
      Coeff.names[2*j] <- paste('Im(',names(X.List)[j],')',sep='')
    }
  }
  names(Model$coefficients) <- Coeff.names
  Coeffs <- as.data.frame(Coeffs,row.names = names(X.List))
  Model$coeffs.complex = Coeffs
  return(Model)
}


## generate random colour --------------
random_color <- function(n){
  colout <- numeric(n)
  for(i in 1:n){
    colorArr = c('1','2','3','4','5','6','7','8','9','A','B','C','D','E','F')
    col <- sample(colorArr,6,replace=T)
    colout[i] <- paste(c('#',col),sep='',collapse='')
  }
  return(colout)
}
