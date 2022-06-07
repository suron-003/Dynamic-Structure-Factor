from re import L
import numpy as np
from numpy import random, linspace, cos, pi
import math
import random
import matplotlib.pyplot as plt
from pyrsistent import b
from scipy.fft import fft, fftfreq
import copy
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits import mplot3d
from plotly import __version__
import pandas as pd
from scipy.optimize import fsolve
import cmath

global l 
global br
global J
global T
global ini_cond
global t_max
global tstep
global mcsteps
global a
global sites
global vec31
global vec32
global vec34
global vec35
global vec36
global vec37
global tarr
global units

l = 6
br = 6
units = l*br
sites = 7*l*br
#anti-ferromagnetic
J = 1
T = 0.1*J
ini_cond = 5
t_max = 200/J
tstep = 0.2/J
tarr = np.arange(0,t_max,tstep)
#site length
a = 1
#vectors connecting different sites in the unit cell
vec31 = np.array([a/2,(a/2)*np.sqrt(3),0])
vec32 = np.array([a,0,0])
vec34 = np.array([a/2,a/(2*np.sqrt(3)),(a/3)*np.sqrt(6)])
vec35 = np.array([a/2,-a/(2*np.sqrt(3)),(2*a*np.sqrt(6))/3])
vec36 = vec35 + np.array([(-a/2),(a/2)*np.sqrt(3),0])
vec37 = vec35 + np.array([(a/2),(a/2)*np.sqrt(3),0])
