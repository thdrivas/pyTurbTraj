import numpy as np
import sympy as sp

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import isotropic1024coarse as info
import pickle
import time

##############################################################################
##############################################################################
##############################################################################
##############################################################################

nu = info['nu']
energy_diss = info['diss']
kolmogorov_time = (nu/energy_diss)**.5
kolmogorov_length = nu**(3/4.)*energy_diss**(-1/4.)

npoints = 10
nparticles = 1000
nsteps = info['time'].shape[0] #total of 124
numcombs =  np.float(nparticles*(nparticles - 1)/2.)

T = info['time'][-1];
subdivisions = 2
t = np.linspace(info['time'][-1], info['time'][0], num = subdivisions*nsteps + 1)
tau = t[0]-t
dt              = t[0] - t[1]
sqrtdt          = abs(dt)**.5

x0 = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
x0[..., 0] = info['lx']*np.random.random(size = (npoints,))[:, None]
x0[..., 1] = info['ynodes'][info['ynodes'].shape[0]//2]
x0[..., 2] = info['lz']*np.random.random(size = (npoints,))[:, None]


trytimes = [1,3,10,30,100,300,1000] #waiting times in case database fails   


pickle.dump(x0, open( "data_isotropic/x0.p", "wb" ) )
#x0 = pickle.load( open( "data_isotropic/x0.p", "rb" ) )

##############################################################################
##############################################################################
##############################################################################
##############################################################################

PrandtlNumbers = np.array([1e1, 1e0, 1e-1])
for m in range(PrandtlNumbers.shape[0]): 
    x    = x0.copy()
    r    = np.zeros(shape = (npoints, nparticles, nparticles),  dtype = np.float32)
    disp = np.zeros(shape = (subdivisions*nsteps + 1, npoints), dtype = np.float32)
    
    Prandtl = np.float(PrandtlNumbers[m])
    kappa = nu/Prandtl
    noiseamplitude  = (2*kappa)**.5
    
    lJHTDB = libJHTDB()
    lJHTDB.initialize()
    for tindex in range(subdivisions*nsteps):
        t0 = time.time()
        print('step {0} of {1} for Pr = {2}'.format(tindex, subdivisions*nsteps, Prandtl))      
        
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
                        r[k,i,j] =    np.sum(np.square(x[k, i, :] - x[k, j, :]))  
        disp[tindex + 1] = np.sum(np.sum(r, axis=2),axis=1)/float(numcombs)
        t1 = time.time()
        time_of_step = t1-t0
        num_steps_left = subdivisions*nsteps - tindex
        time_left = (time_of_step*num_steps_left/60.)/60.
        print('took {0} seconds, about {1} hours remaining'.format(time_of_step,time_left))
    lJHTDB.finalize()

    ####### Dump Data #######
    suffix = 'Traj_{0}_Pr_{1}.p'.format(nparticles, Prandtl)
    pickle.dump(x,    open( "data_isotropic/x"    + suffix, "wb" ) )
    pickle.dump(t,    open( "data_isotropic/t"    + suffix, "wb" ) )
    pickle.dump(disp, open( "data_isotropic/disp" + suffix, "wb" ) )