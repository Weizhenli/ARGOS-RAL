# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 09:56:37 2023

@author: Weizhen Li
"""

import numpy as np
# import matplotlib.pyplot as plt
from scipy.integrate import odeint

def Keller_segel(x, t, init, *args):
    D1, D2, chi, f, k = list(args)[0]
    dx = x[1]-x[0]
    N = len(x)
    kappa = 2*np.pi*np.fft.fftfreq(N, d=dx)
    u0, v0 = init
    dt = t[1]-t[0]
    # solve Keller-­Segel Model
    def rhs_keller1(uv,t,kappa,D1,D2,chi,f,k):
        uv_len = int(len(uv)/2)
        u = uv[:uv_len]
        v = uv[uv_len:]
        
        uhat = np.fft.fft(u)
        vhat = np.fft.fft(v)
        
        d_uhat = (1j)*kappa*uhat
        dd_uhat = -np.power(kappa,2)*uhat
        d_u = np.fft.ifft(d_uhat)
        dd_u = np.fft.ifft(dd_uhat)
        
        d_vhat = (1j)*kappa*vhat
        dd_vhat = -np.power(kappa,2)*vhat
        d_v = np.fft.ifft(d_vhat)
        dd_v = np.fft.ifft(dd_vhat)
        
        du_dt = D1*dd_u - chi*(d_u*d_v + u*dd_v)
        dv_dt = D2*dd_v + f*u - k*v
        
        return np.hstack((du_dt.real,dv_dt.real))

    def rhs_keller2(uv,t,kappa,D1,D2,chi,f,k):
        uv_len = int(len(uv)/2)
        u = uv[uv_len:]
        v = uv[:uv_len]
        
        uhat = np.fft.fft(u)
        vhat = np.fft.fft(v)
        
        d_uhat = (1j)*kappa*uhat
        dd_uhat = -np.power(kappa,2)*uhat
        d_u = np.fft.ifft(d_uhat)
        dd_u = np.fft.ifft(dd_uhat)
        
        d_vhat = (1j)*kappa*vhat
        dd_vhat = -np.power(kappa,2)*vhat
        d_v = np.fft.ifft(d_vhat)
        dd_v = np.fft.ifft(dd_vhat)
        
        du_dt = D1*dd_u - chi*(d_u*d_v + u*dd_v)
        dv_dt = D2*dd_v + f*u - k*v
        
        return np.hstack((dv_dt.real,du_dt.real))
    
    uv = np.hstack((u0,v0))
    out_u = odeint(rhs_keller1, uv, t, args=(kappa,D1,D2,chi,f,k))
    
    vu = np.hstack((v0,u0))
    out_v = odeint(rhs_keller2, vu, t, args=(kappa,D1,D2,chi,f,k))
    
    
    each_len = int(len(uv)/2)
    u_out_u = out_u[:, :each_len]
    v_out_u = out_u[:, each_len:]
    u_out_v = out_v[:, :each_len]
    v_out_v = out_v[:, each_len:]
    
    return u_out_u, v_out_u, u_out_v, v_out_v, dx, dt
    
# L = 10 # Length of domain
# N = 100 # Number of discretization points
# dx = L/N
# x = np.arange(-L/2,L/2,dx) # Define x domain
# kappa = 2*np.pi*np.fft.fftfreq(N, d=dx)
dx = 0.01
x = np.arange(-5, 5, dx)
dt = 0.01
t = np.arange(0,100*dt,dt)

# initial conditions
# u0 = np.sin(2 * np.pi * x)
# v0 = np.cos(2 * np.pi * x)
u0 = 1 + np.zeros(len(x))
# u0 = 1 + np.random.randn(len(x))*1e-2
# u0 = np.exp(-x**2)*0.1
# u0 = np.random.randn(len(x))*1e-3
# v0 = 1 + 1*np.exp(-1*x**2)
v0 = np.exp(-x**2)
# u0 = 1 + np.random.normal(0,0.01,len(x))
# v0 = 0.1*np.exp(-10*x**2)

uv0 = [u0, v0]

coef = [0.1, 1.0, 1.0, 1.0, 1.0]

u_u, v_u, u_v, v_v, dx, dt = Keller_segel(x, t, uv0, coef)



# def rhs_keller(uv,t,kappa,D1,D2,chi,f,k):
#     uv_len = int(len(uv)/2)
#     u = uv[:uv_len]
#     v = uv[uv_len:]
    
#     uhat = np.fft.fft(u)
#     vhat = np.fft.fft(v)
    
#     d_uhat = (1j)*kappa*uhat
#     dd_uhat = -np.power(kappa,2)*uhat
#     d_u = np.fft.ifft(d_uhat)
#     dd_u = np.fft.ifft(dd_uhat)
    
#     d_vhat = (1j)*kappa*vhat
#     dd_vhat = -np.power(kappa,2)*vhat
#     d_v = np.fft.ifft(d_vhat)
#     dd_v = np.fft.ifft(dd_vhat)
    
#     du_dt = D1*dd_u - chi*d_u*d_v - chi*u*dd_v
#     dv_dt = D2*dd_v + f*u - k*v
    
#     return np.hstack((du_dt.real,dv_dt.real))

# constants
# D1 = 1.0 *0.1
# D2 = 1.0
# chi = 1.0 *0.1
# f = 1.0
# k = 1.0

# uv = np.hstack((u0,v0))
# out = odeint(rhs_keller, uv, t, args=(kappa,D1,D2,chi,f,k))

# each_len = int(len(uv)/2)
# u = out[:, :each_len]
# v = out[:, each_len:]

# from matplotlib import pyplot as plt
# plt.imshow(u)
# plt.imshow(v)

# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator

# X, T = np.meshgrid(x, t)
# fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
# surf1 = ax1.plot_surface(T, X, u, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# plt.xlabel('t')
# plt.ylabel('x')

# fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
# surf2 = ax2.plot_surface(T, X, v, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# plt.xlabel('t')
# plt.ylabel('x')
