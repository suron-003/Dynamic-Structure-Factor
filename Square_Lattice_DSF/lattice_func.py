from header import *
def rand_vec():
  a = np.random.standard_normal()
  b = np.random.standard_normal()
  c = np.random.standard_normal()
  l = np.sqrt(a**2+b**2+c**2)
  vec = [a/l,b/l,c/l]

  return vec

#for lattice initialisation
def get_lattice(s):
  s = int(s)
  lat = np.arange(int(s*s*3)).reshape(s,s,3)
  lat = np.double(lat)

  #initialising the lattice with random values
  for i in range (0,s):
    for j in range (0,s):
      lat[i,j] = rand_vec()
  
  return lat

def plot_lattice(lat):
  x = np.array([])
  y = np.array([])
  z = np.array([])
  s = int (len(lat))

  for i in range(0,s):
    for j in range(0,s):
      x = np.append(x,lat[i,j][0])
      y = np.append(y,lat[i,j][1])
      z = np.append(z,lat[i,j][2])

  data = {'x':x,'y':y,'z':z}
  df = pd.DataFrame(data)
  fig = px.scatter_3d(df, x="x",y="y",z="z",width=1000, height=800)
  
  fig.show()

def quiver_plot(lat):
  fig = plt.figure()
  L = int (len(lat))

  x, y = np.meshgrid(np.arange(0, L, 1,dtype=int),
                        np.arange(0, L, 1,dtype=int))
  u = np.array([])
  v = np.array([])
  w = np.array([])

  for i in range(0,L):
    for j in range(0,L):
      u = np.append(u,lat[i,j][0])
      v = np.append(v,lat[i,j][1])
      w = np.append(w,lat[i,j][2])

  u = u.reshape(L,L)
  v = v.reshape(L,L)
  w = w.reshape(L,L)

  plt.quiver(x,y,u,v,w, alpha = 1)
  plt.show()

def save_data(lat,temp):
  s = int (len(lat))
  name = str(s) + "_" + str(temp) + "K"
  filename = "%s.csv" % name
  with open(filename, 'w') as f:
    for i in range(0,s):
      for j in range(0,s):
        string = str (lat[i,j][0]) + "," + str (lat[i,j][1]) + "," + str (lat[i,j][2]) + "\n"
        f.write(string)