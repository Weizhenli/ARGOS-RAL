library(Deriv)
## follow the book Hand book of nonlinear partial differential equations page 858
f <- function(x,t,a1,a2,B1,B2){
  theta1 <- a1*x - a1^3*t
  theta2 <- a2*x - a2^3*t
  A <- ((a1 - a2)/(a1 + a2))^2
  -2*log(1 + B1*exp(theta1) + B2*exp(theta2) + A*B1*B2*exp(theta1+theta2))
}
u_func <- Deriv(Deriv(f,'x'),'x')

n<-512
m=201
kdv.x<-seq(-10,30,length.out = n)
#dt<-0.1
dx<-kdv.x[2]-kdv.x[1]
#t<-seq(dt,m*dt,length.out=m)
kdv.t<-seq(0,20,length.out=m)
#m <- length(t)
dt <- kdv.t[2]-kdv.t[1]

kdv.usol <- matrix(NA,nrow=n,ncol=m)
for(i in 1:n){
  for(j in 1:m){
    kdv.usol[i,j]<-u_func(kdv.x[i],kdv.t[j],0.5,1,1,5)
  }
}

# image(kdv.x,kdv.t,kdv.usol)