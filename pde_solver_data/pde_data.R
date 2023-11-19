setwd('/mnt/d/GitHub/ARGOS-RAL/')
library(reticulate)
source_python('pde_solver_data/ode_int.py')

## advection-diffusion
ad.x = x = seq(-10,10,0.1)
ad.t = t = seq(0,10,0.01)
init = exp(-(x+2)**2)
ad_out = pde_data(c(1,-1),c('u_{xx}','u_{x}'),x,t,init) # advection-diffusion
ad_u <- t(as.matrix(ad_out[[1]]))
dx <- 0.1; dt <- 0.01

ad_list <- list(sol=ad_u,x=ad.x,t=ad.t)

## Cable 
cable.x = x = seq(-4,4,0.1)
cable.t = t = seq(0,5,0.01)
init = exp(-(x)**2)
cable_out = pde_data(c(1,-1),c('u_{xx}','u'),x,t,init) # Cable 
cable_u <- t(as.matrix(cable_out[[1]]))
dx <- 0.1; dt <- 0.01

cable_list <- list(sol=cable_u,x=cable.x,t=cable.t)

# Heat
heat_sol <- function(x,t,k){
  L = max(x)
  u_sol <- matrix(0,nrow=length(x),ncol=length(t))
  for(i in 1:length(x)){
    for(j in 1:length(t)){
      u_sol[i,j] <- 6*sin(pi*x[i]/L)*exp(-k*(pi/L)^2*t[j])
    }
  }
  return(u_sol)
}
dx=0.01;dt=0.01
heat.x = x <- seq(dx,5+dx,dx)
heat.t = t <- seq(0,1.5,dt)
k <- 10

init <- as.vector(heat_sol(x,0,2))
heat_out <- pde_data(list(1),list('u_{xx}'),x,t,init) 
heat_u <- t(as.matrix(heat_out[[1]]))

heat_list <- list(sol=heat_u,x=heat.x,t=heat.t)

# Keller-Segel
source_python('pde_solver_data/Keller_Segel_model.py')
keller.u <- t(u_u)
keller.v <- t(v_u)
keller_sol <- list(u=keller.u,v=keller.v)
keller.x <- x
keller.t <- t

keller_list <- list(sol=keller_sol,x=keller.x,t=keller.t)

# transport
dx = 0.01; dt=0.01
trans.x = x <- seq(-5,1,dx)
trans.t = t <- seq(0,2,dt)
trans_u <- matrix(0,nrow=length(x),ncol=length(t))
for(i in 1:length(x)){
  for(j in 1:length(t)){
    trans_u[i,j] <- exp(-(x[i]+3*t[j])^2)
  }
}

trans_list <- list(sol=trans_u,x=trans.x,t=trans.t)

# qho 
source("pde_solver_data/qho_solver.R")
qho_list <- list(sol=qho.usol,x=xx,t=tt)

# out
save(ad_list,cable_list,heat_list,keller_list,trans_list,qho_list, file='pde_solver_data/pde_data1.RData')
