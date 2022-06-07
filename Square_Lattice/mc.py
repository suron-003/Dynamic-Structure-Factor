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

@jit(nopython=True)
def cell_energy(lat,x,y,vec):
  x1 = int((x+1)%s)
  x2 = int((x-1)%s)
  y1 = int((y+1)%s)
  y2 = int((y-1)%s)
  dxx1 = np.array([0.0,-1.0,0.0])
  dxx2 = np.array([0.0,1.0,0.0])
  dyy2 = np.array([1.0,0.0,0.0])
  dyy1 = np.array([-1.0,0.0,0.0]) 
  b_vec = np.array([0.0,0.0,1.0])
  ehb = J*(np.dot(vec,lat[x1,y]+lat[x2,y]+lat[x,y1]+lat[x,y2]))
  exx1 = -1*d*np.dot(dxx1,np.cross(vec,lat[x1,y]))
  exx2 = -1*d*np.dot(dxx2,np.cross(vec,lat[x2,y]))
  eyy1 = -1*d*np.dot(dyy1,np.cross(vec,lat[x,y1]))
  eyy2 = -1*d*np.dot(dyy2,np.cross(vec,lat[x,y2]))
  edm = (exx1 + exx2 + eyy1 + eyy2)
  eb = -1*b*np.dot(b_vec,vec)

  return ((eb+edm+ehb)) 

@jit(nopython=True)
def mcsweep(lat,T):
  for i in range(0,s):
    for j in range(0,s):
      x = random.randint(0,s-1)
      y = random.randint(0,s-1)
      vec = rand_vec()
      de = cell_energy(lat,x,y,vec) - cell_energy(lat,x,y,lat[x,y])
      
      if (de < 0):
        lat[x,y] = vec
      
      else:
        prob = math.exp((-1*de)/T)
        ran = random.uniform(0,1)
        if (ran < prob):
          lat[x,y] = vec

  return lat 