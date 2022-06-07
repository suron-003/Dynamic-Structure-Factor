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
from parameters import *


def rand_spin():
  a = np.random.standard_normal()
  b = np.random.standard_normal()
  c = np.random.standard_normal()
  l = np.sqrt(a**2+b**2+c**2)
  vec = np.array([a/l,b/l,c/l])

  return vec

def nbrbuild():
#needs to be verified
#constructing the neighbor table
#number of unit cells along the horizontal direction
  #number of unit cells along the vertical direction
  nbr = np.zeros([sites,7])
  
  for i in range(0,sites):
    nbr[i][0] = int (i+1)

  for i in range(0,sites):
    s = int (nbr[i][0])
    if(s%7==1):
      nbr[i][1] = s+1
      nbr[i][2] = s+2
      nbr[i][3] = s+(7*l)+2
      nbr[i][4] = s+(7*l)-6

      if(int(s/(7*l))==(b-1)):
          nbr[i][3] = s-(7*(b-1)*l)+2
          nbr[i][4] = s-(7*(b-1)*l)-6
      if(s%(7*l)==1):
        nbr[i][4] = s+(14*l)-6
        if(int(s/(7*l))==(b-1)):
          nbr[i][4] = (7*l)-5
      nbr[i][5] = s+3
    
    if(s%7==2):
      nbr[i][1] = s-(7*l)+6
      if(s<7*l):
        nbr[i][1] = s + (7*l*(b-1)) + 6

      nbr[i][2] = s+1
      nbr[i][3] = s+8
      if((s+5)%(7*l)==0):
        nbr[i][3] = s - 7*(l-1) + 1
        nbr[i][1] = s - 7*(2*l-1) - 1
        if(s<7*l):
          nbr[i][1] = s + (7*l*(b-2)) + 6
        
      nbr[i][4] = s-1
      nbr[i][5] = s+2

    if(s%7==3):
      nbr[i][1] = s-1
      nbr[i][2] = s - (7*l) - 2
      if(s<(7*l)):
        nbr[i][2] = s + (7*l*(b-1)) - 2
      nbr[i][3] = s-2
      nbr[i][4] = s-8
      if(s%(7*l)==3):
        nbr[i][4] = s + (7*(l-1)) - 1
      nbr[i][5] = s+1

    if(s%7==4):
      nbr[i][1] = s-3
      nbr[i][2] = s-2
      nbr[i][3] = s-1
      nbr[i][4] = s+1
      nbr[i][5] = s+2
      nbr[i][6] = s+3
    
    if(s%7==5):
      nbr[i][1] = s - (7*(l-1)) + 1
      nbr[i][2] = s - (7*l) + 2
      if(s<(7*l)):
        nbr[i][1] = s + (7*l*(b-1)) + 8
        nbr[i][2] = s + (7*l*(b-1)) + 2
      if(s%(7*l)==(7*l)-2):
        nbr[i][1] = s - (7*(2*l-1)) + 1
        if(s<(7*l)):
          nbr[i][1] = 7*b*(l-1) + 6
      nbr[i][3] = s+2
      nbr[i][4] = s+1
      nbr[i][5] = s-1
    
    if(s%7==6):
      nbr[i][1] = s-1
      nbr[i][2] = s-6
      nbr[i][3] = s+1
      nbr[i][4] = s + (7*(l-1)) - 1
      if(s%(7*l)==6):
        nbr[i][2] = s + (7*(l-1)) + 1
        nbr[i][4] = s + (7*(2*l-1)) - 1
      if(int(s/(7*l))==(b-1)):
        nbr[i][4] = s - (7*l*(b-1)) - 8
        if(s%(7*l)==6):
          nbr[i][4] = (7*l) - 2
      nbr[i][5] = s-2
    
    if(s%7==0):
      nbr[i][1] = s+6
      nbr[i][2] = s-2
      nbr[i][3] = s + (7*l) - 2
      nbr[i][4] = s-1
      if(s%(7*l)==0):
        nbr[i][1] = s - (7*(l-1)) - 1
      if(((s/(7*l))>(b-1)) or (int(s/(7*l))==(b))):
        nbr[i][3] = s - (7*(b-1)*l) - 2
      nbr[i][5] = s-3

  #initialising the spins  
  #physics convention for spherical polar coordinates 
  return nbr

def lattice():
    spins = np.zeros([sites,4])
    for i in range(0,sites):
        spins[i][0] = int (i+1)
        spins[i,1:] = rand_spin()

    return spins

