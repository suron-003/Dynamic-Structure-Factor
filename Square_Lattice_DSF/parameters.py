from numpy import tril_indices_from
from header import *
#parameters - for the whole code goes here
global T
global j1
global d
global b
global size

#DM vectors
global dxx1
global dxx2
global dyy1
global dyy2
global gamma

#dynamics
global t_step
global eq_steps
global Th
global Tl
global Tstep

j1 = 1
b = 1*j1
d = 0*j1
size = 10
T = 0.2*j1
t_step = 0.2/j1
eq_steps = size**2
Th = 3*j1
Tl = 0.5*j1
Tstep = 0.5*j1


gamma = 0.65977
dxx1 = d*np.array([0,-1,0]) #*-1 right? [0,1,0] 
dxx2 = d*np.array([0,1,0]) #*-1 right?
dyy2 = d*np.array([1,0,0])
dyy1 = d*np.array([-1,0,0]) 
b_vec = b*np.array([0,0,1])
