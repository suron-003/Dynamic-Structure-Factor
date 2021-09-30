from header import *
from parameters import *


def H_eff(lattice,x,y):
  lat = copy.deepcopy(lattice)  
  s = int (len(lat))  
  x1 = int((x+1)%s)
  x2 = int((x-1)%s)
  y1 = int((y+1)%s)
  y2 = int((y-1)%s)
  heff = (-1*gamma*b_vec) - (0.5*j1*(lat[x1,y] + lat[x2,y] + lat[x,y1] + lat[x,y2])) - (0.5*(np.cross(dxx1,lat[x1,y]) + np.cross(dxx2,lat[x2,y]) + np.cross(dyy1,lat[x,y1]) + np.cross(dyy2,lat[x,y2])))
  return heff

def tstep_euler(lattice): #returns value at t+t_step (t_step = h). Since, for dsf calculations, we need values after every time step.
  lat = copy.deepcopy(lattice)  
  s = int (len(lat))
  for x in range(0,s):
    for y in range(0,s):
      heff = H_eff(lat,x,y)
      dsdt = -1*np.cross(lat[x,y],heff)
      lat[x,y] = lat[x,y] + (t_step*dsdt)
      mag = np.linalg.norm(lat[x,y])
      if(mag>0):
        lat[x,y] = lat[x,y]/mag
  return lat

def tstep_posver(lattice_t,lattice_t1):
  lat_t1 = copy.deepcopy(lattice_t1) 
  lat_t = copy.deepcopy(lattice_t)  
  lat_t2 = copy.deepcopy(lattice_t)
  s = int (len(lat_t))
  for x in range(0,s):
    for y in range(0,s):
      hefft = H_eff(lat_t,x,y,s)
      hefft1 = H_eff(lat_t1,x,y,s)
      lat_t2[x,y] = (2*lat_t[x,y] - lat_t1[x,y]) + ((t_step)*(np.cross(hefft1,lat_t1[x,y]) - np.cross(hefft,lat_t[x,y])))
      mag = np.linalg.norm(lat_t2[x,y])
      if(mag>0):
        lat_t2[x,y] = lat_t2[x,y]/mag
      
  return lat_t2

def time_evol(lat,delt):
    time = np.arange(0,delt,t_step)
    for t in time:
        lat = tstep_euler(lat)