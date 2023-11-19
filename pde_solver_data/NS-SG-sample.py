# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 09:13:24 2021

@author: Administrator
"""

import cv2 as cv
import scipy.signal as ss
import numpy as np
import scipy.io as sio
import itertools
from matplotlib import pyplot as plt
import glob,os,re,time,datetime,sys
import pandas as pd
from sklearn import metrics
#path = r'D:/Users/Administrator/OneDrive - Durham University/Weizhen Li sharing folder/codes/python/'
# path = r'/ddn/data/cfzh32/for_python/TVR'
path = r'/nobackup/cfzh32/for_python/TVR'
sys.path.append(path)
from func.PDE_FIND2_1 import *
import multiprocessing
from functools import partial

def SG_point(u, dx, deg = 3, diff = 1):
    width = len(u)
    index = int((width-1)/2)
    derivatives = (ss.savgol_filter(u, width, deg, diff)/dx**diff)[index]
    return(derivatives)

def mf_point(u, index,kern='Gauss'):
    #weight = np.array([[[0,1,0],[1,2,1],[0,1,0]],[[1,2,1],[2,4,2],[1,2,1]],[[0,1,0],[1,2,1],[0,1,0]]])/28 # Guassian 
    #w = [1,2,1]
    if kern == 'Gauss':
        w = cv.getGaussianKernel(3,0)
        weight = np.outer(np.outer(w, w),w).reshape(3,3,3)
    if kern == 'mean':
        weight = np.repeat(1/27, 27).reshape(3,3,3)
    #weight = weight/np.sum(weight)
    [x,y,t] = index
    u_ = np.sum(weight*u[(x-1):(x+2),(y-1):(y+2),(t-1):(t+2)])
    #weight = np.array([[1,2,1],[2,4,2],[1,2,1]])
    #[x,y,t] = index
    #u_ = np.mean(weight*u[(x-1):(x+2),(y-1):(y+2),t])
    #u_ = cv.GaussianBlur(u[(x-1):(x+2),(y-1):(y+2),(t-1):(t+2)],(3,3,3))
    return(u_)

def mef_point(u,index):
    [x,y,t] = index
    weight = np.repeat(1, 3*3*3).reshape(3,3,3)
    weight[1,1,1] = 0
    u_ = np.median(np.reshape(weight*u[(x-1):(x+2),(y-1):(y+2),(t-1):(t+2)],(3*3*3)))
    return(u_)

def sam_ns_parallel(p,j,d,points_xyt,Wn,dt,dx,dy):
    [x,y,t] = points_xyt[p]
    w_t = SG_point(Wn[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
    w_x = SG_point(Wn[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
    w_y = SG_point(Wn[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
    out = [w_t,w_x,w_y]
    return(out)

def sam_ns(num_xy,num_t,deg=None,noise=0.01,denoise='mf',boundary_max=10,deg_min=3,deg_max=6,cores=2):
    num_xy = int(num_xy); num_t = int(num_t); boundary_max = int(boundary_max)
    deg_min = int(deg_min); deg_max = int(deg_max); cores=int(cores)
    n = 449
    m = 199
    timesteps=151
    dt = 0.2
    dx = 0.02
    dy = 0.02
    
    #path = r'D:/Users/Administrator/OneDrive - Durham University/Weizhen Li sharing folder/codes/python/'
    # path = r'/ddn/data/cfzh32/for_python/TVR/func/'
    path = r'/nobackup/cfzh32/for_python/TVR/func/'
    U1 = np.load(path+'U1.npy')
    V1 = np.load(path+'V1.npy')
    W1 = np.load(path+'W1.npy')
    
    xmin = 100
    xmax = 425
    ymin = 15
    ymax = 185
    
    W = W1[xmin:xmax,ymin:ymax,:]
    U = U1[xmin:xmax,ymin:ymax,:]
    V = V1[xmin:xmax,ymin:ymax,:]
    n,m,steps = W.shape
    np.random.seed(10)
    Wn = W + noise*np.std(W)*np.random.randn(n,m,steps)#- 0.01*np.std(W)*np.random.randn(n,m,steps)
    Un = U + noise*np.std(U)*np.random.randn(n,m,steps)#- 0.01*np.std(U)*np.random.randn(n,m,steps)
    Vn = V + noise*np.std(V)*np.random.randn(n,m,steps)#- 0.01*np.std(V)*np.random.randn(n,m,steps)
    
    if denoise == 'svd':
        Wn = Wn.reshape(n*m,steps)
        Un = Un.reshape(n*m,steps)
        Vn = Vn.reshape(n*m,steps)
        uwn,sigmawn,vwn = np.linalg.svd(Wn, full_matrices=False); vwn = vwn.T
        uun,sigmaun,vun = np.linalg.svd(Un, full_matrices=False); vun = vun.T
        uvn,sigmavn,vvn = np.linalg.svd(Vn, full_matrices=False); vvn = vvn.T
        
        #dim_w = 26;dim_u = 20;dim_v = 20
        dim_w = int(np.array(np.where(np.abs(np.diff(sigmawn))>1e-2))[:,-1])
        dim_u = int(np.array(np.where(np.abs(np.diff(sigmaun))>1e-2))[:,-1])
        dim_v = int(np.array(np.where(np.abs(np.diff(sigmavn))>1e-2))[:,-1])
        
        Ws = uwn[:,0:dim_w].dot(np.diag(sigmawn[0:dim_w]).dot(vwn[:,0:dim_w].T)).reshape(n,m,steps)
        Us = uun[:,0:dim_u].dot(np.diag(sigmaun[0:dim_u]).dot(vun[:,0:dim_u].T)).reshape(n,m,steps)
        Vs = uvn[:,0:dim_v].dot(np.diag(sigmavn[0:dim_v]).dot(vvn[:,0:dim_v].T)).reshape(n,m,steps)
        
        Wn = Wn.reshape(n,m,steps)
        Un = Un.reshape(n,m,steps)
        Vn = Vn.reshape(n,m,steps)
    
    #deg = 5
    #boundary_x = 10
    #np.random.seed(0)
    #num_xy = 300
    #num_t = 5
    #num_points = num_xy * num_t
    if deg == None:
        g_deg = np.arange(deg_min,deg_max+1)
        x_s = np.random.choice(np.arange(boundary_max+1,n-boundary_max-1),num_xy,replace=True)
        y_s = np.random.choice(np.arange(boundary_max+1,m-boundary_max-1),num_xy,replace=True)
        t_s = np.random.choice(np.arange(boundary_max+1,timesteps-boundary_max-1),num_t,replace=False)
        points_xyt = np.vstack(([np.vstack((x_s,y_s,list(itertools.repeat(k,num_xy)))).T for k in t_s]))
        w_gau = np.array([mf_point(Wn,[x,y,t]) for x,y,t in points_xyt]).reshape(num_xy*num_t,1)
        mse_t = np.ones((deg_max,boundary_max))*np.inf
        mse_x = np.ones((deg_max,boundary_max))*np.inf
        mse_y = np.ones((deg_max,boundary_max))*np.inf
        for d in g_deg:
            #bound_min = int((d+7-(d % 2)-1)/2)
            bound_min = int((d+3-(d % 2)-1)/2)
            g_bound = np.arange(bound_min,boundary_max+1)
            for j in g_bound:
                if int(cores)>1:
                    sam_ns_parallel2 = partial(sam_ns_parallel, j=j, d=d, points_xyt=points_xyt, Wn=Wn, dt=dt, dx=dx, dy=dy)
                    #cores = 6
                    pool = multiprocessing.Pool(processes=cores)
                    res1 = pool.map(sam_ns_parallel2, list(range(len(points_xyt))))
                    pool.close()
                    pool.join()
                    res1 = np.array(res1)
                    w_t = res1[:,0]; w_x = res1[:,1]; w_y = res1[:,2]
                else:
                    w_t = np.zeros((len(points_xyt),1))
                    w_x = np.zeros((len(points_xyt),1))
                    w_y = np.zeros((len(points_xyt),1))
                    for p in range(len(points_xyt)):
                        [x,y,t] = points_xyt[p]
                        w_t[p] = SG_point(Wn[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
                        w_x[p] = SG_point(Wn[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
                        w_y[p] = SG_point(Wn[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
                mse_t[d-1,j-1] = metrics.mean_squared_error(w_gau,w_t)
                mse_x[d-1,j-1] = metrics.mean_squared_error(w_gau,w_x)
                mse_y[d-1,j-1] = metrics.mean_squared_error(w_gau,w_y)
        deg_t = np.argmin(pd.DataFrame(mse_t).apply(np.min,axis=1)) + 1
        bound_t = np.argmin(pd.DataFrame(mse_t).apply(np.min,axis=0)) + 1
        deg_x = np.argmin(pd.DataFrame(mse_x).apply(np.min,axis=1)) + 1
        bound_x = np.argmin(pd.DataFrame(mse_x).apply(np.min,axis=0)) + 1
        deg_y = np.argmin(pd.DataFrame(mse_y).apply(np.min,axis=1)) + 1
        bound_y = np.argmin(pd.DataFrame(mse_y).apply(np.min,axis=0)) + 1
        
    else:
        deg[1] = int(deg[1]); deg[0] = int(deg[0]); 
        boundary_xy = deg[1]
        boundary_t = deg[0]
        x_s = np.random.choice(np.arange(boundary_xy+1,n-boundary_xy-1),num_xy,replace=True)
        y_s = np.random.choice(np.arange(boundary_xy+1,n-boundary_xy-1),num_xy,replace=True)
        t_s = np.random.choice(np.arange(boundary_t+1,timesteps-boundary_t-1),num_t,replace=False)
        points_xyt = np.vstack(([np.vstack((x_s,y_s,list(itertools.repeat(i,num_xy)))).T for i in t_s]))
        deg_t_u = deg[0]
        deg_x = deg_y = deg[1]
        bound_t = boundary_t
        bound_x = bound_y = boundary_xy
        
    w = np.zeros((len(points_xyt),1))
    u = np.zeros((len(points_xyt),1))
    v = np.zeros((len(points_xyt),1))
    wt = np.zeros((len(points_xyt),1))
    wx = np.zeros((len(points_xyt),1))
    wy = np.zeros((len(points_xyt),1))
    wxx = np.zeros((len(points_xyt),1))
    wxy = np.zeros((len(points_xyt),1))
    wyy = np.zeros((len(points_xyt),1))
    
    for p in range(len(points_xyt)):
        #p=1
        [x,y,t] = points_xyt[p]
        if denoise == 'mf':
            w[p] = mf_point(Wn,[x,y,t])#Ws[x,y,t]
            u[p] = mf_point(Un,[x,y,t])#Us[x,y,t]
            v[p] = mf_point(Vn,[x,y,t])#Vs[x,y,t]
            wt[p] = SG_point([mf_point(Wn,[x,y,t0]) for t0 in range((t-bound_t),(t+bound_t+1))],dt,deg=deg_t,diff=1) 
            wx[p] = SG_point([mf_point(Wn,[x0,y,t]) for x0 in range((x-bound_x),(x+bound_x+1))],dx,deg=deg_x,diff=1)
            wy[p] = SG_point([mf_point(Wn,[x,y0,t]) for y0 in range((y-bound_y),(y+bound_y+1))],dy,deg=deg_y,diff=1)
            wxx[p] = SG_point([mf_point(Wn,[x0,y,t]) for x0 in range((x-bound_x),(x+bound_x+1))],dx,deg=deg_x,diff=2) 
            wyy[p] = SG_point([mf_point(Wn,[x,y0,t]) for y0 in range((y-bound_y),(y+bound_y+1))],dy,deg=deg_y,diff=2) 
            x_diff_yp = SG_point([mf_point(Wn,[x0,y+1,t]) for x0 in range((x-bound_x),(x+bound_x+1))],dx,deg=deg_x,diff=1) 
            x_diff_ym = SG_point([mf_point(Wn,[x0,y-1,t]) for x0 in range((x-bound_x),(x+bound_x+1))],dx,deg=deg_x,diff=1) 
            wxy[p] = (x_diff_yp-x_diff_ym)/(2*dy)
        if denoise == 'svd':
            w[p] = Ws[x,y,t]
            u[p] = Us[x,y,t]
            v[p] = Vs[x,y,t]
            wt[p] = SG_point(Ws[x,y,(t-boundary_t):(t+boundary_t+1)],dt,deg=deg,diff=1) 
            wx[p] = SG_point(Ws[(x-bound_t):(x+bound_t+1),y,t],dx,deg=deg,diff=1) 
            wy[p] = SG_point(Ws[x,(y-bound_y):(y+bound_y+1),t],dy,deg=deg,diff=1)
            wxx[p] = SG_point(Ws[(x-bound_x):(x+bound_x+1),y,t],dx,deg=deg,diff=2) 
            wyy[p] = SG_point(Ws[x,(y-bound_y):(y+bound_y+1),t],dy,deg=deg,diff=2) 
            x_diff_yp = SG_point(Ws[(x-bound_y):(x+bound_y+1),y+1,t],dx,deg=deg,diff=1) 
            x_diff_ym = SG_point(Ws[(x-bound_y):(x+bound_y+1),y-1,t],dx,deg=deg,diff=1) 
            wxy[p] = (x_diff_yp-x_diff_ym)/(2*dy)
        
    X_data = np.hstack([w,u,v])
    X_ders = np.hstack([np.ones((len(points_xyt),1)), wx, wy, wxx, wxy, wyy])
    X_ders_descr = ['','w_{x}', 'w_{y}','w_{xx}','w_{xy}','w_{yy}']
    X, description = build_Theta(X_data, X_ders, X_ders_descr, 2, data_description = ['w','u','v'])
    df = np.real(np.hstack((wt,X)))
    df = pd.DataFrame(df)
    df.columns = ['wt']+description
    
    return(df)




def sam_rd_parallel(p,j,d,points_xyt,Un,Vn,dt,dx,dy):
    [x,y,t] = points_xyt[p]
    u_t = SG_point(Un[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
    u_x = SG_point(Un[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
    u_y = SG_point(Un[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
    v_t = SG_point(Vn[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
    v_x = SG_point(Vn[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
    v_y = SG_point(Vn[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
    #out = np.vstack((u_t,u_x,u_y,v_t,v_x,v_y)).T
    out = [u_t,u_x,u_y,v_t,v_x,v_y]
    return(out)


def sam_rd(num_xy,num_t,deg=None,noise=0.005,boundary_max=10,deg_min=3,deg_max=6,cores=2):
    #data = sio.loadmat('F:/reaction_diffusion_big.mat')
    num_xy = int(num_xy); num_t = int(num_t); boundary_max = int(boundary_max)
    deg_min = int(deg_min); deg_max = int(deg_max); cores = int(cores)
    # data = sio.loadmat('/ddn/home/cfzh32/for_matlab/reaction_diffusion_big.mat')
    data = sio.loadmat('/nobackup/cfzh32/for_matlab/reaction_diffusion_big.mat')
    t = data['t'][:,0]
    x = data['x'][0,:]
    y = data['y'][0,:]
    U = data['u']
    V = data['v']
    
    np.random.seed(10)
    n = len(x) # also the length of y
    steps = len(t)
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    dt = t[2]-t[1]
    Un = U + noise*np.std(U)*np.random.randn(n,n,steps)
    Vn = V + noise*np.std(V)*np.random.randn(n,n,steps)
    
    #deg = 5
    #boundary = deg*2
    #boundary_x = 10
    #np.random.seed(0)
    #num_xy = 500
    #num_t = 5*2
    
    if deg == None:
        g_deg = np.arange(deg_min,deg_max+1)
        x_s = np.random.choice(np.arange(boundary_max+1,n-boundary_max-1),num_xy,replace=True)
        y_s = np.random.choice(np.arange(boundary_max+1,n-boundary_max-1),num_xy,replace=True)
        t_s = np.random.choice(np.arange(boundary_max+1,steps-boundary_max-1),num_t,replace=False)
        points_xyt = np.vstack(([np.vstack((x_s,y_s,list(itertools.repeat(k,num_xy)))).T for k in t_s]))
        u_gau = np.array([mf_point(Un,[x,y,t]) for x,y,t in points_xyt]).reshape(num_xy*num_t,1)
        v_gau = np.array([mf_point(Vn,[x,y,t]) for x,y,t in points_xyt]).reshape(num_xy*num_t,1)
        mse_t_u = np.ones((deg_max,boundary_max))*np.inf
        mse_x_u = np.ones((deg_max,boundary_max))*np.inf
        mse_y_u = np.ones((deg_max,boundary_max))*np.inf
        mse_t_v = np.ones((deg_max,boundary_max))*np.inf
        mse_x_v = np.ones((deg_max,boundary_max))*np.inf
        mse_y_v = np.ones((deg_max,boundary_max))*np.inf
        for d in g_deg:
            #bound_min = int((d+7-(d % 2)-1)/2)
            bound_min = int((d+3-(d % 2)-1)/2)
            g_bound = np.arange(bound_min,boundary_max+1)
            for j in g_bound:
                if int(cores)>1:
                    sam_rd_parallel2 = partial(sam_rd_parallel, j=j, d=d, points_xyt=points_xyt, Un=Un, Vn=Vn, dt=dt, dx=dx, dy=dy)
                    #cores = 6
                    pool = multiprocessing.Pool(processes=cores)
                    res1 = pool.map(sam_rd_parallel2, list(range(len(points_xyt))))
                    pool.close()
                    pool.join()
                    #[u_t,u_x,u_y,v_t,v_x,v_y] = res1
                    res1 = np.array(res1)
                    u_t = res1[:,0]; u_x = res1[:,1]; u_y = res1[:,2]
                    v_t = res1[:,3]; v_x = res1[:,4]; v_y = res1[:,5]
                else:
                    u_t = np.zeros((len(points_xyt),1))
                    u_x = np.zeros((len(points_xyt),1))
                    u_y = np.zeros((len(points_xyt),1))
                    v_t = np.zeros((len(points_xyt),1))
                    v_x = np.zeros((len(points_xyt),1))
                    v_y = np.zeros((len(points_xyt),1))
                    for p in range(len(points_xyt)):
                        [x,y,t] = points_xyt[p]
                        u_t[p] = SG_point(Un[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
                        u_x[p] = SG_point(Un[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
                        u_y[p] = SG_point(Un[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
                        v_t[p] = SG_point(Vn[x,y,(t-j):(t+j+1)],dt,deg=d,diff=0)
                        v_x[p] = SG_point(Vn[(x-j):(x+j+1),y,t],dx,deg=d,diff=0)
                        v_y[p] = SG_point(Vn[x,(y-j):(y+j+1),t],dy,deg=d,diff=0)
                mse_t_u[d-1,j-1] = metrics.mean_squared_error(u_gau,u_t)
                mse_x_u[d-1,j-1] = metrics.mean_squared_error(u_gau,u_x)
                mse_y_u[d-1,j-1] = metrics.mean_squared_error(u_gau,u_y)
                mse_t_v[d-1,j-1] = metrics.mean_squared_error(v_gau,v_t)
                mse_x_v[d-1,j-1] = metrics.mean_squared_error(v_gau,v_x)
                mse_y_v[d-1,j-1] = metrics.mean_squared_error(v_gau,v_y)
        
        deg_t_u = np.argmin(pd.DataFrame(mse_t_u).apply(np.min,axis=1)) + 1
        bound_t_u = np.argmin(pd.DataFrame(mse_t_u).apply(np.min,axis=0)) + 1
        deg_x_u = np.argmin(pd.DataFrame(mse_x_u).apply(np.min,axis=1)) + 1
        bound_x_u = np.argmin(pd.DataFrame(mse_x_u).apply(np.min,axis=0)) + 1
        deg_y_u = np.argmin(pd.DataFrame(mse_y_u).apply(np.min,axis=1)) + 1
        bound_y_u = np.argmin(pd.DataFrame(mse_y_u).apply(np.min,axis=0)) + 1
        
        deg_t_v = np.argmin(pd.DataFrame(mse_t_v).apply(np.min,axis=1)) + 1
        bound_t_v = np.argmin(pd.DataFrame(mse_t_v).apply(np.min,axis=0)) + 1
        deg_x_v = np.argmin(pd.DataFrame(mse_x_v).apply(np.min,axis=1)) + 1
        bound_x_v = np.argmin(pd.DataFrame(mse_x_v).apply(np.min,axis=0)) + 1
        deg_y_v = np.argmin(pd.DataFrame(mse_y_v).apply(np.min,axis=1)) + 1
        bound_y_v = np.argmin(pd.DataFrame(mse_y_v).apply(np.min,axis=0)) + 1
    else:
        deg[1] = int(deg[1]); deg[0] = int(deg[0]); 
        boundary_xy = deg[1]
        boundary_t = deg[0]
        x_s = np.random.choice(np.arange(boundary_xy+1,n-boundary_xy-1),num_xy,replace=True)
        y_s = np.random.choice(np.arange(boundary_xy+1,n-boundary_xy-1),num_xy,replace=True)
        t_s = np.random.choice(np.arange(boundary_t+1,steps-boundary_t-1),num_t,replace=False)
        points_xyt = np.vstack(([np.vstack((x_s,y_s,list(itertools.repeat(i,num_xy)))).T for i in t_s]))
        deg_t_u = deg_t_v = deg[0]
        deg_x_u = deg_y_u = deg_x_v = deg_y_v = deg[1]
        bound_t_u = bound_t_v = boundary_t
        bound_x_u = bound_y_u = bound_x_v = bound_y_v = boundary_xy
    
    u = np.zeros((len(points_xyt),1),dtype=complex)
    v = np.zeros((len(points_xyt),1),dtype=complex)
    ut = np.zeros((len(points_xyt),1),dtype=complex)
    vt = np.zeros((len(points_xyt),1),dtype=complex)
    ux = np.zeros((len(points_xyt),1),dtype=complex)
    uy = np.zeros((len(points_xyt),1),dtype=complex)
    uxx = np.zeros((len(points_xyt),1),dtype=complex)
    uxy = np.zeros((len(points_xyt),1),dtype=complex)
    uyy = np.zeros((len(points_xyt),1),dtype=complex)
    vx = np.zeros((len(points_xyt),1),dtype=complex)
    vy = np.zeros((len(points_xyt),1),dtype=complex)
    vxx = np.zeros((len(points_xyt),1),dtype=complex)
    vxy = np.zeros((len(points_xyt),1),dtype=complex)
    vyy = np.zeros((len(points_xyt),1),dtype=complex)
    
    for p in range(len(points_xyt)):
        [x,y,t] = points_xyt[p]
        u[p] = mf_point(Un,[x,y,t])
        v[p] = mf_point(Vn,[x,y,t])
        # time derivatives
        ut[p] = SG_point([mf_point(Un,[x,y,t0]) for t0 in range((t-bound_t_u),(t+bound_t_u+1))],dt,deg=deg_t_u,diff=1)
        vt[p] = SG_point([mf_point(Vn,[x,y,t0]) for t0 in range((t-bound_t_v),(t+bound_t_v+1))],dt,deg=deg_t_v,diff=1) 
            
        # spatial derivatives
        ux_diff_yp = SG_point([mf_point(Un,[x0,y+1,t]) for x0 in range((x-bound_x_u),(x+bound_x_u+1))],dx,deg=deg_x_u,diff=1)
        ux_diff_ym = SG_point([mf_point(Un,[x0,y-1,t]) for x0 in range((x-bound_x_u),(x+bound_x_u+1))],dx,deg=deg_x_u,diff=1)
        vx_diff_yp = SG_point([mf_point(Vn,[x0,y+1,t]) for x0 in range((x-bound_x_v),(x+bound_x_v+1))],dx,deg=deg_x_v,diff=1)
        vx_diff_ym = SG_point([mf_point(Vn,[x0,y-1,t]) for x0 in range((x-bound_x_v),(x+bound_x_v+1))],dx,deg=deg_x_v,diff=1)
        
        ux[p] = SG_point([mf_point(Un,[x0,y,t]) for x0 in range((x-bound_x_u),(x+bound_x_u+1))],dx,deg=deg_x_u,diff=1)
        uy[p] = SG_point([mf_point(Un,[x,y0,t]) for y0 in range((y-bound_y_u),(y+bound_y_u+1))],dy,deg=deg_y_u,diff=1) 
        uxx[p] = SG_point([mf_point(Un,[x0,y,t]) for x0 in range((x-bound_x_u),(x+bound_x_u+1))],dx,deg=deg_x_u,diff=2)
        uxy[p] = (ux_diff_yp-ux_diff_ym)/(2*dy)
        uyy[p] = SG_point([mf_point(Un,[x,y0,t]) for y0 in range((y-bound_y_u),(y+bound_y_u+1))],dy,deg=deg_y_u,diff=2)
        
        vx[p] = SG_point([mf_point(Vn,[x0,y,t]) for x0 in range((x-bound_x_v),(x+bound_x_v+1))],dx,deg=deg_x_v,diff=1)
        vy[p] = SG_point([mf_point(Vn,[x,y0,t]) for y0 in range((y-bound_y_v),(y+bound_y_v+1))],dy,deg=deg_y_v,diff=1) 
        vxx[p] = SG_point([mf_point(Vn,[x0,y,t]) for x0 in range((x-bound_x_v),(x+bound_x_v+1))],dx,deg=deg_x_v,diff=2)
        vxy[p] = (vx_diff_yp-vx_diff_ym)/(2*dy)
        vyy[p] = SG_point([mf_point(Vn,[x,y0,t]) for y0 in range((y-bound_y_v),(y+bound_y_v+1))],dy,deg=deg_y_v,diff=2)
    X_data = np.hstack([u,v])
    X_ders = np.hstack([np.ones((len(points_xyt),1)), ux, uy, uxx, uxy, uyy, vx, vy, vxx, vxy, vyy])
    X_ders_descr = ['','u_{x}', 'u_{y}','u_{xx}','u_{xy}','u_{yy}','v_{x}', 'v_{y}','v_{xx}','v_{xy}','v_{yy}']
    X, description = build_Theta(X_data, X_ders, X_ders_descr, 3, data_description = ['u','v'])
    df = np.real(np.hstack((ut,vt,X)))
    df = pd.DataFrame(df)
    df.columns = ['ut','vt']+description
    
    return(df)
