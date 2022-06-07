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


@jit(nopython=True)
def ene(lat):
  e = 0
  for x in range(0,s):
    for y in range(0,s):
      x1 = int((x+1)%s)
      x2 = int((x-1)%s)
      y1 = int((y+1)%s)
      y2 = int((y-1)%s)
      ej = (J*(np.dot(lat[x,y],lat[x1,y]+lat[x2,y]+lat[x,y1]+lat[x,y2]))) 
      exx1 = -1*d*np.dot(dxx1,np.cross(lat[x,y],lat[x1,y]))
      exx2 = -1*d*np.dot(dxx2,np.cross(lat[x,y],lat[x2,y]))
      eyy1 = -1*d*np.dot(dyy1,np.cross(lat[x,y],lat[x,y1]))
      eyy2 = -1*d*np.dot(dyy2,np.cross(lat[x,y],lat[x,y2]))
      edm = (exx1 + exx2 + eyy1 + eyy2)
      eb = -1*b*np.dot(b_vec,lat[x,y])
      e = e + (ej/2 + edm/2 + eb) 
  
  return e

@jit(nopython=True)
def H_eff(lat,x,y):
  x1 = int((x+1)%s)
  x2 = int((x-1)%s)
  y1 = int((y+1)%s)
  y2 = int((y-1)%s)
  dxx1 = np.array([0.0,-1.0,0.0])
  dxx2 = np.array([0.0,1.0,0.0])
  dyy2 = np.array([1.0,0.0,0.0])
  dyy1 = np.array([-1.0,0.0,0.0]) 
  b_vec = np.array([0.0,0.0,1.0])
  heff = (-1*b*gamma*b_vec) + (0.5*J*(lat[x1,y] + lat[x2,y] + lat[x,y1] + lat[x,y2])) - d*(0.5*(np.cross(dxx1,lat[x1,y]) + np.cross(dxx2,lat[x2,y]) + np.cross(dyy1,lat[x,y1]) + np.cross(dyy2,lat[x,y2])))
  return heff