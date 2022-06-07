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

global J
global s
global stepsize
global ini_cond
global T
global d
global b
global dxx1
global dxx2
global dyy1
global dyy2
global gamma
global tarr
global t_max
global mcsteps

#ferromagnetic
J = -1
s = 4
b = 0.0
d = 0.0
stepsize = 0.01/np.abs(J)
t_max = 10/np.abs(J)
tarr = np.arange(0,t_max,stepsize)
ini_cond = 2
Tc = 0.69*np.abs(J)
T = 0.4*np.abs(J)
mcsteps = s**2
gamma = 0.65977
dxx1 = np.array([0,-1.0,0])
dxx2 = np.array([0,1.0,0])
dyy2 = np.array([1.0,0,0])
dyy1 = np.array([-1.0,0,0]) 
b_vec = np.array([0,0,1.0])