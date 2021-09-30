from numpy import ravel_multi_index
from mc import *
from parameters import *
from header import *
from dynamics import *
from lattice_func import *

def sqt(s,ini_cond,q_array,delt):#takes it size, q vector and max time and returns S(q,t) for all q and all t
    s = int (s)
    #time = np.arange(0,delt,t_step)
    time = np.arange(0,delt + t_step,t_step)
    sqt_mat = (np.zeros([int (len(q_array)),int (len(time))]))
    for i in range(0,ini_cond+1):
        print(i)
        lat = get_lattice(s)
        #for jjj in range(0,200):
        #    lat = mc_sweep(lat,Tl)
        lat = annl(lat,Th,Tl,Tstep,eq_steps)
        sq0 = (np.zeros([int(len(q_array)),3]))
        for q in range(0,len(q_array)):
            qx = q_array[q][0]
            qy = q_array[q][1]
            for x in range(0,s):
                for y in range(0,s):
                    sq0[q] = sq0[q] + ((cmath.exp(1j*((qx*(x-(0.5*s)) + qy*(y-(0.5*s))))))*lat[x,y])
        
        sq0 = sq0/s
        for t in range(0,len(time)):
            sqat = np.zeros([len(q_array),3])
            batq = np.zeros(len(q_array))
            for q in range(0,len(q_array)):
                qx = q_array[q][0]
                qy = q_array[q][1]
                for x in range(0,s):
                    for y in range(0,s):
                        sqat[q] = sqat[q] + ((cmath.exp(1j*((qx*(x-(0.5*s)) + qy*(y-(0.5*s))))))*lat[x,y])
                        #((cmath.exp(1j*((qx*x + qy*y))))*lat[x,y])
                sqat[q] = sqat[q]/s
                batq[q] = np.abs(np.dot(sqat[q],np.conjugate(sq0[q])))
            
            #for ppp in range(0,len(sqat)):
            #    sqt_mat[ppp,t] = sqt_mat[ppp,t] + sqat[ppp]
            sqt_mat[:,t] = sqt_mat[:,t] + batq
            #print(lat)
            lat = tstep_euler(lat)
            #print('Hello')
            #print(lat)
    sqt_mat = sqt_mat/ini_cond
    for ii in range(0,len(q_array)):
        sqt_mat[ii,:] = sqt_mat[ii,:]/sqt_mat[ii][0]
    
    return sqt_mat

def ft(st):
    sw = np.array([])
    for k in range(0,len(st)):
        s0 = 0
        for i in range(0,len(st)):
            s0 = s0 + st[i]*complex(np.cos(k*i*((2*np.pi)/len(st))), -1*np.sin(k*i*((2*np.pi)/len(st))))
        sw = np.append(sw,s0)
    return sw

                 
def ssc(s,ini_cond,q_array,delt): #spin-spin correlator
    time = np.arange(0, delt + t_step, t_step)
    cmatt = np.zeros([int (len(time)),s**2,s**2])
    for i in range(0,ini_cond):
        lat = get_lattice(s)
        lat_t = copy.deepcopy(lat)
        ctemp = np.zeros([s**2,s**2])
        for t in range(0,len(time)):
            for p in range(0,s**2):
                for q in range(0,s**2):
                    rip = int (p/s) 
                    rjp = p % s
                    riq = int (q/s)
                    rjq = q % s
                    ctemp[p,q] = np.dot(lat[rip,rjp],lat[riq,rjq])

            cmatt[t] = cmatt[t] + ctemp
    
    cmatt = cmatt/ini_cond

    return cmatt
        
def dsf_analytic(lattice,q_array,delt):
    lat = copy.deepcopy(lattice)
    s = int (len(lat))
    time = np.arange(0,delt + t_step,t_step)
    sqt_mat = (np.zeros([int (len(q_array)),int (len(time))]))
    sq0 = (np.zeros([int(len(q_array)),3]))
    for q in range(0,len(q_array)):
        qx = q_array[q][0]
        qy = q_array[q][1]
        for x in range(0,s):
            for y in range(0,s):
                sq0[q] = sq0[q] + ((cmath.exp(1j*((qx*(x-(0.5*s)) + qy*(y-(0.5*s))))))*lat[x,y])
    
    sq0 = sq0/s
    for t in range(0,len(time)):
        sqat = np.zeros([len(q_array),3])
        batq = np.zeros(len(q_array))
        for q in range(0,len(q_array)):
            qx = q_array[q][0]
            qy = q_array[q][1]
            for x in range(0,s):
                for y in range(0,s):
                    sqat[q] = sqat[q] + ((cmath.exp(1j*((qx*(x-(0.5*s)) + qy*(y-(0.5*s))))))*lat[x,y])
                    #((cmath.exp(1j*((qx*x + qy*y))))*lat[x,y])
            sqat[q] = sqat[q]/s
            batq[q] = np.abs(np.dot(sqat[q],np.conjugate(sq0[q])))
        
        sqt_mat[:,t] = sqt_mat[:,t] + batq
        lat = tstep_euler(lat)

    for ii in range(0,len(q_array)):
        sqt_mat[ii,:] = sqt_mat[ii,:]/sqt_mat[ii][0]
    return sqt_mat






