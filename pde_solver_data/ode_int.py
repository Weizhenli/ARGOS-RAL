# %%
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import re
from scipy.fftpack import diff as psdiff

def p_diff_term(u, kappa, diff):
    diff = int(diff)
    uhat = np.fft.fft(u)
    d_uhat = np.power((1j)*kappa,diff)*uhat
    d_u = np.fft.ifft(d_uhat)
    return(d_u)

def p_diff_term(u, dx, N, diff):
    diff = int(diff)
    L = dx*(N)
    d_u = psdiff(u, period=L, order=diff)
    return(d_u)


def poly_diff_order(string):
    string_new = re.sub('u_','',string)
    diff_order = len(re.findall('x', string_new))
    if len(re.findall('\\^', string_new)) == 0:
        poly_order = len(re.findall('u', string_new))
    else:
        poly_order = int(re.search('\d', string_new).group())
    return(int(poly_order), int(diff_order))

def pde_data(coefs, term_names, x, t, init, method='RK45'):
    # nu = 0.1 # Diffusion constant
    # L = 20     # Length of domain
    # N = 1000   # Number of discretization points
    # dx = L/N
    # x = np.arange(-L/2,L/2,dx) # Define x domain
    dx = x[1]-x[0]
    N = len(x)

    # Define discrete wavenumbers
    kappa = 2*np.pi*np.fft.fftfreq(N, d=dx)

    # Initial condition
    # u0 = 1/np.cosh(x)
    # u0 = np.exp(-(x+2)**2)
    u0 = init
    dt = t[1]-t[0]
    # Simulate PDE in spatial domain
    # dt = 0.01
    # t = np.arange(0,100*dt,dt)

    # def pde(u,t,kappa,coefs,term_names):
    #     uhat = np.fft.fft(u)
    #     du_dt = 0
    #     for i in range(len(coefs)):
    #         poly, diff = poly_diff_order(term_names[i])
    #         if poly == 1 and diff == 0:
    #             poly = 0; diff = 0
    #         # du_dt += coefs[i]*np.power(u,poly)*p_diff_term(u, kappa, diff)
    #         du_dt += coefs[i]*np.power(u,poly)*p_diff_term(u, dx, N, diff)
    #     return du_dt.real

    # u = odeint(pde,u0,t,args=(kappa,coefs,term_names))

    def pde(t,u,kappa,coefs,term_names):
        uhat = np.fft.fft(u)
        du_dt = 0
        for i in range(len(coefs)):
            poly, diff = poly_diff_order(term_names[i])
            if poly == 1 and diff == 0:
                poly = 0; diff = 0
            # du_dt += coefs[i]*np.power(u,poly)*p_diff_term(u, kappa, diff)
            du_dt += coefs[i]*np.power(u,poly)*p_diff_term(u, dx, N, diff)
        return du_dt.real

    u = solve_ivp(pde, [t[0],t[-1]], u0, args=(kappa,coefs,term_names), t_eval=t, method=method).y.T

    return(u, dx, dt)

# %% 
# if __name__ == '__main__':
#     import matplotlib.pyplot as plt
#     from matplotlib import cm
#     ## Burgers 
#     x = np.arange(-10,10,0.1)
#     t = np.arange(0,10,0.01)
#     xx, tt = np.meshgrid(x, t)
#     init = np.exp(-(x+2)**2)
#     u, dx, dt = pde_data([-1.0589,0.158214],['uu_{x}','u_{xx}'],x,t,init) # burgers
#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#     ax.plot_surface(xx, tt, u.reshape((len(t),len(x))), rstride=1, cstride=1,cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
#     ## kdv
#     def kdv_exact(x, c):
#         """Profile of the exact solution to the KdV for a single soliton on the real line."""
#         u = 0.5*c*np.cosh(0.5*np.sqrt(c)*x)**(-2)
#         return u
#     dx = 0.15
#     x = np.arange(-20,20,dx)
#     t = np.arange(0,10,0.01)
#     xx, tt = np.meshgrid(x, t)
#     init = kdv_exact(x+5,1) + kdv_exact(x-8, 0.5)
#     # init = kdv_exact(x-0.33*L, 0.75) + kdv_exact(x-0.65*L, 0.4)
#     u, dx, dt = pde_data([-1,-6],['u_{xxx}','uu_{x}'],x,t,init) # kdv
#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#     ax.plot_surface(xx, tt, u.reshape((len(t),len(x))), rstride=1, cstride=1,cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
#     ## heat equation 
#     u, dx, dt = pde_data([-1.0589],['u_{xx}'],x,t,init)
    
#     import os, sys
#     ## add path to the folder containing the ode_int.py file
#     sys.path.append('C:/Users/cfzh32/Documents/GitHub/Weizhen-Li/predict')
#     from deepxde_funcs import *
#     dt = t[1]-t[0]
#     dx = x[2]-x[1]

#     def init_fun(x):
#         return(np.exp(-(x[:,0:1]+2)**2))
#     u_pred, f_2 = predict_uxt([-1.0589],['u_{xx}'],x,t,init_fun,3000)

#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#     ax.plot_surface(xx, tt, u_pred.reshape((len(t),len(x))), rstride=1, cstride=1,cmap=cm.coolwarm, linewidth=0, antialiased=False)

# %%
