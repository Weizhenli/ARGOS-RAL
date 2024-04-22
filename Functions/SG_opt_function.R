library(signal)
library(Deriv)
library(Metrics)
library(purrr)
library(pracma)
library(reticulate)
cv <- import("cv2")

my_SG_filter <- function(x, p, ws, dx, d = 0){
  p <- round(p,0)
  ws <- round(ws,0)
  if(p < 2){
    p <- 2
  }
  # ws <- as.integer(param[2])
  if(ws > length(x)-1){
    ws = length(x)-1
  }
  if(ws <= p){
    ws = p + 1
  }
  if(!(ws %% 2)){  # then make odd
    ws <-  ws + 1
  }
  x_smooth <- sgolayfilt(x, p=p, n=ws, m=d, ts=dx)
  return(x_smooth)
}


object_fun <- function(param, dx, y_noise, y_true) {
  p <- param[1]
  ws <- param[2]
  y_smooth <- my_SG_filter(y_noise, p, ws, dx)
  obj <- rmse(y_smooth, y_true)
  # return(sqrt(sum((y_smooth-y_true)^2)))
  return(obj)
}

my_SG_filter_matrix <- function(u, p, ws, dx, d = 0){
  p <- round(p,0)
  ws <- round(ws,0)
  if(p < 2){
    p <- 2
  }
  # ws <- as.integer(param[2])
  if(ws > nrow(u)-1){
    ws = nrow(u)-1
  }
  if(ws <= p){
    ws = p + 1
  }
  if(!(ws %% 2)){  # then make odd
    ws <-  ws + 1
  }
  u_smooth <- apply(u, 2, sgolayfilt, p = p, n = ws, m = d, ts = dx)
  return(u_smooth)
}
object_fun_mat <- function(param, dx, y_noise, y_true) {
  p <- param[1]
  ws <- param[2]
  y_smooth <- my_SG_filter_matrix(y_noise, p, ws, dx)
  obj <- rmse(y_smooth, y_true)
  # return(sqrt(sum((y_smooth-y_true)^2)))
  return(obj)
}

correct_params <- function(params, params_low, params_high){
  sapply(seq_along(params), function(i){
    param <- round(params[i],0)
    param = max(param, params_low[i])
    param = min(param, params_high[i])
    param = round(param,0)
    return(param)
  })
}
opt_func <- function(param, fun, y_noise, y_true, dx=dx, lower, upper){
  result <-
    tryCatch(
      result <- optim(param, fn = fun, y_noise = y_noise, y_true = y_true, dx=dx)[1:2],
      error = function(e) {
        par <- correct_params(param, lower, upper)
        value <- fun(par, dx=dx, y_noise=y_noise, y_true=y_true)
        return(list(par=par, value=value))
      }
    )
  result$par <- round(result$par)
  return(result)
}

opt_SG_param <- function(y_noise, y_true, dx, object_fun){
  params <- list()
  # orders = c(2, 4, 8, 10, 15)
  # window_sizes = c(3, 10, 30, 50, 90, 130, 200, 300)
  orders = c(4, 8, 12)
  window_sizes = c(5, 10, 30, 50, 100)
  for(order in orders){
    for(ws in window_sizes){
      params <- c(params, list(c(p=order, ws=ws)))
    }
  }
  
  lower = c(4,5); upper = c(Inf, Inf)
  results <- params |> map(\(x) opt_func(x, fun=object_fun, y_noise=y_noise, y_true=y_true, dx=dx, lower=lower, upper=upper))
  opt_params <- lapply(results, function(x) x$par)
  opt_values <- sapply(results, function(x) x$value)
  opt_index <- which.min(opt_values)
  par = opt_params[[opt_index]]
  if(par[1] < 4){
    par[1] <- 4
  }
  if(!is.null(nrow(y_noise))){
    if(par[2] > length(nrow(y_noise))-1){
      par[2] <- length(nrow(y_noise))-1
    }
  }else{
    if(par[2] > length(y_noise)-1){
      par[2] <- length(y_noise)-1
    }
  }
  if(par[2] <= par[1]){
    par[2] <- par[1] + 1
  }
  if(!(par[2] %% 2)){  # then make odd
    par[2] <-  par[2] + 1
  }
  return(list(par = par, metric = min(opt_values)))
}

## Gaussian blur
GB <- function(x,kernel){
  kernel <- kernel / sum(kernel)
  x2 <- c(x[1],x,x[length(x)])
  convolve(x2, kernel, type="f") 
}

ASG_build_library_opt <- function(u, dx, dt, sg_para=NULL, noise=0.01, diff_o=3, poly_o=3,
                              ref_data='Gblur', ...){
  u <- as.matrix(u)
  if(noise > 0){
    un <- u + rnorm(nrow(u)*ncol(u), 0, noise*sd(u))
  }else{
    un <- u
    #un_gblur <- u
  }
  if(ref_data == 'Gblur'){
    un_gblur <- cv$GaussianBlur(un, as.integer(c(3,3)), 0)
  }else{
    un_gblur <- un
  }
  n=nrow(u);m=ncol(u)
  if(is.null(sg_para)){
    opt_SG_x <- opt_SG_param(y_noise=u, y_true=un_gblur, dx=dx, object_fun=object_fun_mat)
    opt_SG_t <- opt_SG_param(y_noise=t(u), y_true=t(un_gblur), dx=dt, object_fun=object_fun_mat)
    
    deg_x <- opt_SG_x$par[1]
    width_x <- opt_SG_x$par[2]
    
    deg_t <- opt_SG_t$par[1]
    width_t <- opt_SG_t$par[2]
  }else{
    deg_x <- sg_para$deg_x # index of row
    width_x <- sg_para$width_x # index of column
    deg_t <- sg_para$deg_t # index of row
    width_t <- sg_para$width_t # index of column
  }
  
  linear_sys <- build_linear_system(un_gblur, dt, dx, D=diff_o, P=poly_o, time_diff = 'poly', space_diff = 'poly',
                                    deg_x=deg_x, deg_t=deg_t,width_x=width_x, width_t=width_t)
  return(linear_sys)
}

