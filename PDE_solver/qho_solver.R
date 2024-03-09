library(pracma)
L <- 15 # number of domain
N <- 2^9 # number of discretization points
dx <- L/N
xx <- seq(-L/2, L/2-dx, dx)

# Define discrete wavenumbers
kappa = (2*pi/L)*c((-N/2):(N/2-1))
kappa = fftshift(kappa)    # Re-order fft wavenumbers

# Initial condition
u0 = exp(-((xx-1)/2)^2)
# u0 = 2*sech(xx-1)
u = u0

# Simulate PDE in spatial domain
dt = 0.025
tt = seq(0,10,dt)

qho.usol = matrix(0, nrow = length(tt), ncol=length(xx))
qho.usol[1,] = as.complex(u)

for(j in seq_along(tt)[-1]){
  # u1 = u*exp(-1i*(0.5*xx^2)*dt)
  # u2 = ifft(fft(Conj(u1))*exp(1i*kappa^2*dt/2))
  # u = Conj(u2)
  u1 = ifft(fft(u)*exp(1i*kappa^2*dt/2))
  u2 = u1*Conj(exp(-1i*0.5*(xx^2)*dt))
  u = u2
  qho.usol[j,] = u
}

# image(abs(qho.usol))
# image(Re(qho.usol))
# image(Im(qho.usol))
