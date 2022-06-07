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

def energy(nbr,spins):
  e = 0
  sp = copy.deepcopy(spins)
  for i in range(0,len(nbr)):
    spi = np.array(sp[i,1:])
    nbrs = np.array(nbr[i,1:])
    for n in nbrs:
      n = int(n)
      if(n==0):
        continue
      nspi = np.array(sp[n-1,1:])
      e = e + (0.5)*np.dot(nspi,spi)
  return e