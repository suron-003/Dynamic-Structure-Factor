from header import *
from calc import *

def rand_vec():
  a = np.random.standard_normal()
  b = np.random.standard_normal()
  c = np.random.standard_normal()
  l = np.sqrt(a**2+b**2+c**2)
  vec = [a/l,b/l,c/l]

  return vec

def mc_sweep(lattice,temp):
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  for i in range(0,s):
    for j in range(0,s):
      x = random.randint(0,s-1)
      y = random.randint(0,s-1)
      vec = rand_vec()
      vec1 = copy.deepcopy(lat[x,y])
      de = cell_energy(lat,x,y,vec) - cell_energy(lat,x,y,vec1)
      
      if (de <= 0):
        lat[x,y] = copy.deepcopy(vec)
      
      else:
        prob = math.exp((-1*de)/temp)
        ran = random.uniform(0,1)
        if (ran <= prob):
          lat[x,y] = copy.deepcopy(vec)

  return lat

def annl(lattice,Th,Tl,Tstep,eq_steps):
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  temp_array = np.arange(Th,Tl-Tstep,-1*Tstep)
  for T in temp_array:
    for i in range(0,eq_steps):
      lat = mc_sweep(lat,T)
  return lat
    