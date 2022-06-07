import numpy as np
from numpy import random, linspace, cos, pi
import math
import random
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.fft import rfft, rfftfreq
import copy
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits import mplot3d
from plotly import __version__
import pandas as pd
from scipy.optimize import fsolve
import cmath
from numba import jit
from numpy import linalg as LA
from scipy.linalg import expm, norm
from scipy.integrate import odeint
import time
import numba
from parameters import *
from lattice import *

import plotly.offline as pyo
import plotly.graph_objs as go
from plotly.offline import iplot
import plotly.figure_factory as ff
import plotly.express as px
from mc import *
from obs import *

#dynamics
@jit(nopython=True)
def ode(r,t):
  lat = r.reshape(s,s,3)
  der = np.zeros((s,s,3))
  #der = copy.deepcopy(lat)
  for i in range(0,s):
    for j in range(0,s):
      heff = H_eff(lat,i,j)
      der[i,j] = -1*np.cross(lat[i,j],heff)
  der = der.flatten()

  return der

#correlation functions
@jit(nopython=True)
def sq(lat,q):
  sqt = np.array([0.0+0.0j,0.0+0.0j,0.0+0.0j])
  for x in range(0,s):
    for y in range(0,s):
      r = np.array([x-(s/2),y-(s/2)]) 
      sqt = sqt + ((lat[x,y]*cmath.exp(1j*np.dot(q,r))))
  sqt = (1/s)*sqt

  return sqt

@jit
def dcfps(qmat):
  dqt = np.zeros((len(qmat),len(time)),dtype=complex) 
  sq0 = np.zeros((len(qmat),3),dtype=complex)
  for m in range(0,ini_cond):
    print(m)
    lat = get_lattice(s)
    for j in range(0,mcsteps):
      lat = mcsweep(lat,T)
    for j in range(0,len(qmat)):
      q = np.array(qmat[j,:])
      sq0[j,:] = sq(lat,q)
      dqt[j][0] = dqt[j][0] + np.dot(sq0[j,:],np.conjugate(sq0[j,:]))
    rt = odeint(ode,lat.flatten(),time)
    #latt = copy.deepcopy(lat)
    for t in range(1,len(time)):
      latt = rt[t,:]
      #latt = euler(latt)
      latt = latt.reshape(s,s,3)
      for j in range(0,len(qmat)):
        q = np.array(qmat[j,:])
        dqt[j][t] = dqt[j][t] + np.dot(sq(latt,q),np.conjugate(sq0[j,:]))
  dqt = dqt/ini_cond
  dqt = np.abs(dqt)
  for j in range(0,len(qmat)):
    mag = dqt[j][0]*1
    dqt[j,:] = dqt[j,:]/mag
  
  return dqt


@jit
def vanhove():
  g = np.zeros((len(time),(2*s)-1,(2*s)-1))
  latarr = np.array([0,0])
  for i in range(0,s):
    for j in range(0,s):
      latarr = np.vstack([latarr,np.array([i,j])])
  latarr = latarr[1:,:]
  for i in range(0,ini_cond):
    print(i)
    lat0 = get_lattice(s)
    for p in range(0,mcsteps):
      lat0 = mcsweep(lat0,T)
    rt = odeint(ode,lat0.flatten(),time)
    for t in range(0,len(time)):
      lat1 = rt[t,:]
      lat1 = lat1.reshape(s,s,3)
      for x in range(-s+1,s):
        for y in range(-s+1,s):
          for k in range(0,len(latarr)):
            lark = latarr[k]
            larp = lark + np.array([x,y])
            larp = larp%s
            g[t,x+s-1,y+s-1] = g[t,x+s-1,y+s-1] + (np.dot(lat1[larp[0],larp[1]],lat0[lark[0],lark[1]])/(s**2))
  g = g/ini_cond

  return g

@jit
def dcfvh(qmat):
  g = vanhove()
  dcfmat = np.zeros((len(qmat),len(time)),dtype=complex)
  for i in range(0,len(qmat)):
    qarr = np.array(qmat[i,:])
    for t in range(0,len(time)):
      for x in range(-s+1,s):
        for y in range(-s+1,s):
          dcfmat[i][t] = dcfmat[i][t] + ((g[t,x+s-1,y+s-1])*(cmath.exp(-1j*np.dot(qarr,np.array([x,y])))))

  dcfmat = np.abs(dcfmat)
  
  for j in range(0,len(qmat)):
    mag = dcfmat[j][0]*1
    dcfmat[j,:] = dcfmat[j,:]/mag
  
  return dcfmat

#static structure factor
def ssf(lat):
  kx = np.arange((-1*np.pi), (np.pi) + 0.3, (2*np.pi)/s)
  ky = np.arange((-1*np.pi), (np.pi) + 0.3, (2*np.pi)/s)
  a = len(kx)
  sk=np.array([])
  sk = sk.astype(complex)

  for k1 in kx:
    for k2 in ky:
      sum = 0
      q = np.array([k1,k2])
      for x1 in range(0,s):
        for y1 in range(0,s):
          for x2 in range(0,s):
            for y2 in range(0,s):
              r = np.array([x2-x1,y2-y1])
              sum = sum + (np.dot(lat[x1,y1],lat[x2,y2])*cmath.exp(1j*np.dot(q,r))) 
      sum = sum/(2*s) 
      sum = sum * np.conjugate(sum)
      sk = np.append(sk, np.abs(sum))

  sk = np.abs(sk)
  sk = np.reshape(sk, (a,a))
  fig = go.Figure(data=[go.Surface(z=sk, x=kx, y=ky)])
  fig.update_layout(title='Static Structure Factor', autosize=False,
                    width=500, height=500,
                    margin=dict(l=65, r=50, b=65, t=90))
  fig.show()


def dsfsz(qmat):
    g = vanhove()
    dsfmat = np.zeros((len(qmat),int(len(time)/2)+1),dtype=complex)
    for i in range(0,len(qmat)):
        qarr = np.array(qmat[i,:])
        for x in range(-s+1,s):
            for y in range(-s+1,s):
                darr = np.array([])
                for t in range(0,len(time)):
                    darr = np.append(darr,((g[t,x+s-1,y+s-1])*(cmath.exp(-1j*np.dot(qarr,np.array([x,y]))))))
        dsfmat[i,:] = dsfmat[i,:] + np.rfft(darr)
    dsfmat = dsfmat/(2*np.pi*s**2) 
    dsfmat = np.abs(dsfmat)
    return dsfmat