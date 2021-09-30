from header import *
from lattice_func import *
from parameters import *
from mc import *
from dsf import *
from plot import *


#initialising and annealing
#lat = get_lattice(size)
lat=np.loadtxt('config_one_skyrmion_radius_4.dat')
print(lat[:,2])
print(lat[:,3])
#for x in range(0,size):
#    for y in range(0,size):
#        lat[x,y] = np.array([1,0,0])
lattice = get_lattice(8)
for x in range(0,7):
    for y in range(0,7):
        lattice[x,y] = np.array([np.sin(lat[:,2][x])*np.cos(lat[:,3][y]),

q_arr = np.array([[0,0],[0.2*np.pi,0.2*np.pi],[0.4*np.pi,0.4*np.pi]])
ini_cond = 10
delt = 5/j1
sqt_mat = dsf_analytic(lat,q_arr,delt)
print(sqt_mat)
sqt1 = sqt_mat[0,:]
sqt2 = sqt_mat[1,:]
sqt3 = sqt_mat[2,:]
time = np.arange(0,delt + t_step,t_step)

plt.plot(time,sqt1)
plt.plot(time,sqt2)
plt.plot(time,sqt3)
plt.xlabel('Time')
plt.ylabel('S(t)')
plt.title('S(t) vs t')
plt.legend(["q = 0","q = 0.2*pi","q = 0.4*pi"])
plt.savefig("S_vs_t.png") 
plt.show()

'''
for i in range(0,200):
    lat = mc_sweep(lat,T)

#plot_lattice(lat)
#quiver_plot(lat)

m = calcmag(lat)
print(m)
'''