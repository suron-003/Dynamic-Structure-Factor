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

#lattice functions
@jit(nopython=True)
def rand_vec():
  a = np.random.standard_normal()
  b = np.random.standard_normal()
  c = np.random.standard_normal()
  l = np.sqrt(a**2+b**2+c**2)
  vec = np.array([a/l,b/l,c/l])

  return vec

#@numba.cfunc("array(float32, 3d, A)(int64)")
@jit
def get_lattice(s):
  lat = np.arange(int(s*s*3)).reshape(s,s,3)
  lat = np.double(lat)

  for i in range (0,s):
    for j in range (0,s):
      lat[i,j] = rand_vec()
  
  return lat

def save_data(lat,temp):
  name = str(s) + "_" + str(temp) + "K"
  filename = "%s.csv" % name
  with open(filename, 'w') as f:
    for i in range(0,s):
      for j in range(0,s):
        string = str (lat[i,j][0]) + "," + str (lat[i,j][1]) + "," + str (lat[i,j][2]) + "\n"
        f.write(string)