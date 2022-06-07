import numpy as np
from numpy import random, linspace, cos, pi
import math
import random
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
import copy
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits import mplot3d
from plotly import __version__
import pandas as pd
import numba
from scipy.optimize import fsolve
import cmath
from scipy.integrate import odeint
from numba import jit
from lattice import *
from parameters import *
from mc import *

def ode(r,t):
  lat = r.reshape(sites,4)
  der = np.zeros((sites,4))
  #der = copy.deepcopy(lat)
  nbr = nbrbuild()
  for i in range(0,len(nbr)):
    nbrs = np.array(nbr[i,1:])
    spi = np.array(lat[i,1:])
    dsdt = 0
    for n in nbrs:
      n = int(n)
      if(n==0):
        continue
      nspi = np.array(lat[n-1,1:])
      dsdt = dsdt - (J*np.cross(spi,nspi))
    
    der[i][0] = 0
    der[i,1:] = dsdt

  der = der.flatten()
  return der

def sq(q,spins):
  sqa = np.zeros((3),dtype=np.complex_)
  i = 3
  while(i<(7*l*b)):
    j = i-1
    row = int (i/(7*l))
    col = int ((i%(7*l))/7)
    #base position of the unit cell (position of 3rd site)
    bpos = 2*col*vec32 + 2*row*vec31
    s3 = np.array(spins[j,1:])
    s1 = np.array(spins[j-2,1:])
    s2 = np.array(spins[j-1,1:])
    s4 = np.array(spins[j+1,1:])
    s5 = np.array(spins[j+2,1:])
    s6 = np.array(spins[j+3,1:])
    s7 = np.array(spins[j+4,1:])

    sqa = sqa + np.array((s3*cmath.exp(1j*np.dot(q,bpos)))) + np.array((s1*cmath.exp(1j*np.dot(q,vec31+bpos)))) + np.array((s2*cmath.exp(1j*np.dot(q,vec32+bpos)))) 
    + np.array((s4*cmath.exp(1j*np.dot(q,vec34+bpos)))) + np.array((s5*cmath.exp(1j*np.dot(q,vec35+bpos)))) + np.array((s6*cmath.exp(1j*np.dot(q,vec36+bpos)))) 
    + np.array((s7*cmath.exp(1j*np.dot(q,vec37+bpos)))) 

    #sq = sq + np.array((spins[j,1:]*cmath.exp(1j*np.dot(q,bpos)))) + np.array((spins[j-2,1:]*cmath.exp(1j*np.dot(q,vec31+bpos)))) + np.array((spins[j-1,1:]*cmath.exp(1j*np.dot(q,vec32+bpos)))) 
    #+ np.array((spins[j+1,1:]*cmath.exp(1j*np.dot(q,vec34+bpos)))) + np.array((spins[j+2,1:]*cmath.exp(1j*np.dot(q,vec35+bpos)))) + np.array((spins[j+3,1:]*cmath.exp(1j*np.dot(q,vec36+bpos)))) 
    #+ np.array((spins[j+4,1:]*cmath.exp(1j*np.dot(q,vec37+bpos))))    

    i = i + 7
  sqa = sqa/np.sqrt(7*l*b)

  return sqa

def dcf(qmat,T):
  dcfmat = np.zeros((len(qmat),int(tarr/tstep)),dtype=np.complex_)
  sqa = np.zeros((len(qmat),3),dtype=np.complex_)
  for i in range(0,ini_cond):
    nbr = nbrbuild() 
    spins = lattice()
    #equilibriation
    for st in range(0,sites):
      spins = mc_sweep(nbr,spins,T)
    for qq in range(0,len(qmat)):
      qarr = qmat[qq,:]
      qarr = np.array(qarr)
      sqa[qq] = (sq(qarr,spins))
      dcfmat[qq][0] = dcfmat[qq][0] + (np.dot(sqa[qq],np.conjugate(sqa[qq])))
    
    rt = odeint(ode,spins.flatten(),tarr)
    for tt in range(1,int(tarr/tstep)):
      spins = rt[tt,:]
      spins = spins.reshape((sites,3))
      for qq in range(0,len(qmat)):
        qarr = qmat[qq,:]
        sqt = (sq(qarr,spins))
        dcfmat[qq][tt] = dcfmat[qq][tt] + (np.dot(sqt,np.conjugate(sqa[qq])))
  
  dcfmat = dcfmat/ini_cond
  dcfmat = np.abs(dcfmat)
  for i in range(0,len(dcfmat)):
    mag = dcfmat[i][0]*1
    dcfmat[i,:] = dcfmat[i,:]/(mag)
  
  return dcfmat
