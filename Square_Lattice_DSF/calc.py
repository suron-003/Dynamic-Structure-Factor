from header import *
from parameters import *

def calcmag(lattice):
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  m = 0
  mag = np.array([0,0,0])
  for i in range(0,s):
    for j in range(0,s):
      mag = mag + lat[i,j]
  print(mag)
  m = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
  m = m/s**2
  return (m)

#calculates cell energy
def cell_energy(lattice,x,y,vec):
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  x1 = int((x+1)%s)
  x2 = int((x-1)%s)
  y1 = int((y+1)%s)
  y2 = int((y-1)%s)
  ehb =  (-1*j1*(np.dot(vec,lat[x1,y]+lat[x2,y]+lat[x,y1]+lat[x,y2])))
  exx1 = -1*np.dot(dxx1,np.cross(vec,lat[x1,y]))
  exx2 = -1*np.dot(dxx2,np.cross(vec,lat[x2,y]))
  eyy1 = -1*np.dot(dyy1,np.cross(vec,lat[x,y1]))
  eyy2 = -1*np.dot(dyy2,np.cross(vec,lat[x,y2]))
  edm = (exx1 + exx2 + eyy1 + eyy2)
  eb = -1*np.dot(b_vec,vec)

  return ((eb+edm+ehb))  

#calculate energy of a lattice
def calcenergy_total(lattice):
  e0 = 0
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  for x in range(0,s):
    for y in range(0,s):
      vec = copy.deepcopy(lat[x,y])
      x1 = int((x+1)%s)
      x2 = int((x-1)%s)
      y1 = int((y+1)%s)
      y2 = int((y-1)%s)
      ehb =  (-1*j1*(np.dot(vec,lat[x1,y]+lat[x2,y]+lat[x,y1]+lat[x,y2])))
      exx1 = -1*np.dot(dxx1,np.cross(vec,lat[x1,y]))
      exx2 = -1*np.dot(dxx2,np.cross(vec,lat[x2,y]))
      eyy1 = -1*np.dot(dyy1,np.cross(vec,lat[x,y1]))
      eyy2 = -1*np.dot(dyy2,np.cross(vec,lat[x,y2]))
      edm = (exx1 + exx2 + eyy1 + eyy2)
      eb = -1*np.dot(b_vec,vec)
      e0 = e0 + (ehb/2 + edm/2 + eb) 
  
  return (e0/s**2)

def static_structure_factor(lattice):
  lat = copy.deepcopy(lattice)
  s = int (len(lat))
  kx = np.arange((-1*np.pi), (np.pi) + 0.3, (2*np.pi)/s)
  ky = np.arange((-1*np.pi), (np.pi) + 0.3, (2*np.pi)/s)
  a = len(kx)
  sk=np.array([])

  for k1 in kx:
    for k2 in ky:
      sum = 0
      for x1 in range(0,s):
        for y1 in range(0,s):
          for x2 in range(0,s):
            for y2 in range(0,s):
              #sum = sum + (np.dot(lat[x1,y1],lat[x2,y2])- complex(lat[x1,y1][2]*lat[x2,y2][2]))*complex((np.cos(k1*(x2-x1) + k2*(y2-y1))),(np.sin(k2*(y2-y1) + k1*(x2-x1)))) #S_xy
              sum = sum + np.dot(lat[x1,y1],lat[x2,y2]) * complex((np.cos(k1*(x2-x1) + (k2*(y2-y1)))),(np.sin(k1*(x2-x1) + k2*(y2-y1))))
      #sum = sum/(2*s) #factor of 2 difference max
      sum = sum/2
      sum = sum * np.conjugate(sum)
      sk = np.append(sk, np.abs(sum))

  sk = np.reshape(sk, (a,a))
  fig = go.Figure(data=[go.Surface(z=sk, x=kx, y=ky)])
  fig.update_layout(title='Static Structure Factor', autosize=False,
                    width=500, height=500,
                    margin=dict(l=65, r=50, b=65, t=90))
  fig.show()