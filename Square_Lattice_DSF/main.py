from header import *
from lattice_func import *
from parameters import *
from mc import *
from dsf import *
from plot import *


#initialising and annealing
lat = get_lattice(size)


q_arr = np.array([[0,0],[0.2*np.pi,0.2*np.pi],[0.4*np.pi,0.4*np.pi]])
ini_cond = 1
delt = 50/j1
sqt_mat = sqt(size,ini_cond,q_arr,delt)
sqt1 = sqt_mat[0,:]
sqt2 = sqt_mat[1,:]
sqt3 = sqt_mat[2,:]
time = np.arange(0,delt + t_step,t_step)

'''
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
sw1 = ft(sqt1)
sw2 = ft(sqt2)
sw3 = ft(sqt3)


w = np.linspace((0*np.pi)/(len(sqt1)*t_step),(2*np.pi)/(len(sqt1)*t_step),num = len(sqt1), endpoint = True)
sw_1 = copy.deepcopy(sw1[0: int (len(sw1)/2)])
sw1 = np.append(sw_1[::-1],sw_1)
sw1 = np.append(np.array([0]),sw1)
sw_2 = copy.deepcopy(sw2[0: int (len(sw1)/2)])
sw2 = np.append(sw_2[::-1],sw_2)
sw2 = np.append(np.array([0]),sw2)
sw_3 = copy.deepcopy(sw3[0: int (len(sw1)/2)])
sw3 = np.append(sw_3[::-1],sw_3)
sw3 = np.append(np.array([0]),sw3)
plt.plot(w,sw1)
plt.plot(w,sw2)
plt.plot(w,sw3)
plt.xlabel('omega')
plt.ylabel('S(w)')
plt.title('S(w) vs w')
plt.legend(["q = 0","q = 0.2*pi","q = 0.4*pi"])
#plt.savefig("S_vs_t.png") 
plt.show()