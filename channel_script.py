import numpy as np
import sympy as sp

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import channel as info
import pickle
import time

##############################################################################
##############################################################################
##############################################################################
##############################################################################

nu = 5e-5

Lx = info['lx']
Ly = info['ly']
Lz = info['lz']

Top = Ly/2
Bottom = -Ly/2
eps = Ly/100


npoints = 2
nparticles = 100
numcombs   = np.float(nparticles*(nparticles-1)/2)
subdivisions = 2

T = 25.9935
database_dt = 0.0065
nsteps = int(T/database_dt) #info['time'].shape[0]
t = np.linspace(T, info['time'][0], num = subdivisions*nsteps+1)
tau = t[0]-t
dt              = t[0] - t[1]
sqrtdt          = abs(dt)**.5

x0         = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
x0[..., 0] = Lx*np.random.random(size = (npoints,))[:, None]
x0[..., 1] = 0*info['ynodes'][info['ynodes'].shape[0]//10]
x0[..., 2] = Lz*np.random.random(size = (npoints,))[:, None]

trytimes = [1,3,10,30,100,300,1000] #waiting times in case database fails   

pickle.dump(x0, open( "data_channel/x0.p", "wb" ) )
#x0 = pickle.load( open( "data_channel/x0.p", "rb" ) )

##############################################################################
##############################################################################
##############################################################################
##############################################################################

PrandtlNumbers = np.array([1e1, 1e0, 1e-1])
for m in range(PrandtlNumbers.shape[0]):
    x = x0.copy() 
    LB         = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
    LT         = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
    HT         = (info['time'][-1]+1)**2*np.ones(shape = (npoints, nparticles, 1), dtype = np.float32)
    r          = np.zeros(shape = (npoints, nparticles,  nparticles), dtype = np.float32)
    disp       = np.zeros(shape = (subdivisions*nsteps + 1, npoints), dtype = np.float32)
    
    Prandtl = np.float(PrandtlNumbers[m])
    kappa = nu/Prandtl
    noiseamplitude  = (2*kappa)**.5
    
    lJHTDB = libJHTDB()
    lJHTDB.initialize()
    for tindex in range(subdivisions*nsteps):
        t0 = time.time()
        print('step {0} of {1} for Pr = {2}'.format(tindex,subdivisions*nsteps, Prandtl))
        
        for tryT in trytimes:  
            try:
                u = lJHTDB.getData(
                            t[tindex],
                            x,
                            sinterp = interpolation_code['M2Q8'],
                            tinterp = interpolation_code['NoTInt'],
                            data_set = info['name'],
                            getFunction = 'getVelocity')
                break
            except Exception as e:
                print e
                time.sleep(tryT)
                    
        dW = np.random.randn(*x.shape)*sqrtdt
        dX = -u*dt + noiseamplitude*dW
        x += dX
        for k in range(npoints):
            for i in range(nparticles):
                for j in range(nparticles):
                    if i<j:
                        r[k,i,j] =  np.sum(np.square(x[k, i, :] - x[k, j, :]))           
        disp[tindex + 1] = np.sum(np.sum(r, axis=2),axis=1)/float(numcombs)

        cond1 = x[..., 1] < Bottom + eps
        cond2 = x[..., 1] > Top    - eps
        HT[np.logical_or(cond1, cond2)] = np.minimum(HT[np.logical_or(cond1, cond2)], tau[tindex])
        c1indices = np.where(cond1)   
        if c1indices[0].size > 0:
            V  =  -2*dt*np.log(np.random.random(x.shape))[c1indices]
            Y  = (-dX[c1indices] + np.sqrt(2*kappa*V + dX[c1indices]**2))/2.0
            dL =  Y - x[c1indices] + Bottom
            dL[np.where(dL < 0)] = 0
            x[c1indices]  += dL
            LB[c1indices] += dL/kappa
        c2indices = np.where(cond2) 
        if c2indices[0].size > 0:
            V  =  -2*dt*np.log(np.random.random(x.shape))[c2indices]
            Y  = (-dX[c2indices] + np.sqrt(2*kappa*V + dX[c2indices]**2))/2.0
            dL =  Y + x[c2indices] - Top
            dL[np.where(dL < 0)] = 0
            x[c2indices]  -= dL
            LT[c2indices] += dL/kappa
        x[..., 1] = np.clip(x[..., 1], Bottom, Top)
        t1 = time.time()
        time_of_step = t1-t0
        num_steps_left = subdivisions*nsteps - tindex
        time_left = (time_of_step*num_steps_left/60.)/60.
        print('took {0} seconds, about {1} hours remaining'.format(time_of_step,time_left))
    lJHTDB.finalize() 
    
    ####### Dump Data #######
    suffix = 'middle_Traj_{0}_Pr_{1}.p'.format(nparticles, Prandtl)
    pickle.dump(x, open( "data_channel/x"          + suffix, "wb" ) )
    pickle.dump(t, open( "data_channel/t"           + suffix, "wb" ) )
    pickle.dump(LT, open( "data_channel/LT"         + suffix, "wb" ) )
    pickle.dump(LB, open( "data_channel/LB"         + suffix, "wb" ) )
    pickle.dump(HT, open( "data_channel/HT"         + suffix, "wb" ) )
    pickle.dump(disp, open( "data_channel/disp"     + suffix, "wb" ) )