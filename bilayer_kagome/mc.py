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
from numba import jit
from lattice import *
from parameters import *

#plain heisenberg model H = j Si.Sj
def mc_sweep(nbr,spins,temp):
  sp = copy.deepcopy(spins)
  for i in range(0,len(nbr)):
    rsite = random.randint(1,len(nbr))
    spi = np.array(sp[rsite-1,1:])
    rspin = rand_spin()
    #calculating change in energy
    e1 = 0
    e2 = 0
    nbrs = np.array(nbr[rsite-1,1:])
    for n in nbrs:
      n = int(n)
      if(n==0):
        continue
      nspi = np.array(sp[n-1,1:])
      e1 = e1 + J*np.dot(nspi,spi)
      e2 = e2 + J*np.dot(nspi,rspin)
      #e1 = e1 + J*((np.sin(spi[0])*np.cos(spi[1])*np.sin(nspi[0])*np.cos(nspi[1]))+(np.sin(spi[0])*np.sin(spi[1])*np.sin(nspi[0])*np.sin(nspi[1])) + (np.cos(spi[0])*np.cos(nspi[0])))
      #e2 = e2 + J*((np.sin(rspin[0])*np.cos(rspin[1])*np.sin(nspi[0])*np.cos(nspi[1]))+(np.sin(rspin[0])*np.sin(rspin[1])*np.sin(nspi[0])*np.sin(nspi[1])) + (np.cos(rspin[0])*np.cos(nspi[0])))
    de = e2 - e1

    if (de < 0):
      sp[rsite-1,1:] = rspin
      #sp[rsite-1][1],sp[rsite-1][2] = rspin[0],rspin[1]
      
    else:
      prob = math.exp((-1*de)/temp)
      ran = random.uniform(0,1)
      if (ran < prob):
        sp[rsite-1,1:] = rspin
        #sp[rsite-1][1],sp[rsite-1][2] = rspin[0],rspin[1]
      else:
        continue
  
  return sp