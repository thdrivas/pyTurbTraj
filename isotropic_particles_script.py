import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import isotropic1024coarse as info

import pickle
import time

npoints = 10
nparticles = 100
nsteps = info['time'].shape[0] #total of 124

#x0 = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
#x0[..., 0] = info['lx']*np.random.random(size = (npoints,))[:, None]
#x0[..., 1] = info['ynodes'][info['ynodes'].shape[0]//2]
#x0[..., 2] = info['lz']*np.random.random(size = (npoints,))[:, None]

nu = info['nu']
energy_diss = info['diss']

kolmogorov_time = (nu/energy_diss)**.5
kolmogorov_length = nu**(3/4.)*energy_diss**(-1/4.)

T = info['time'][-1];
subdivisions = 2
t = np.linspace(info['time'][-1], info['time'][0], num = subdivisions*nsteps + 1)
tau = t[0]-t
dt              = t[0] - t[1]
sqrtdt          = abs(dt)**.5

trytimes = [1,3,10,30,100,300,1000]

x0 = pickle.load( open( "data_isotropic/x0.p", "rb" ) )
x0 = x0[:,0:nparticles,:]

PrandtlNumbers = np.array([1e1, 1e0, 1e-1])
for m in range(PrandtlNumbers.shape[0]): 
    xfull  = np.zeros(shape = (subdivisions*nsteps+1, npoints, nparticles, 3), dtype = np.float32)
    xfull[0] = x0.copy()
    
    Prandtl = np.float(PrandtlNumbers[m])
    kappa = nu/Prandtl
    noiseamplitude  = (2*kappa)**.5
    
    lJHTDB = libJHTDB()
    lJHTDB.initialize()
    for tindex in range(subdivisions*nsteps):
        print('step {0} of {1} for Pr = {2}'.format(tindex, subdivisions*nsteps, Prandtl))
        
        for tryT in trytimes:  
            try:
                u = lJHTDB.getData(
                            t[tindex],
                            xfull[tindex],
                            sinterp = interpolation_code['M2Q8'],
                            tinterp = interpolation_code['NoTInt'],
                            data_set = info['name'],
                            getFunction = 'getVelocity')
                break
            except Exception as e:
                print e
                time.sleep(tryT)      
                          
        dW = np.random.randn(*xfull.shape[1:])*sqrtdt
        dX = -u*dt + noiseamplitude*dW
        xfull[tindex + 1] = xfull[tindex] + dX
    lJHTDB.finalize()

    ####### Dump Data #######
    suffix = 'Traj_{0}_Pr_{1}.p'.format(nparticles, Prandtl)
    pickle.dump(xfull,    open( "data_isotropic/xfull"    + suffix, "wb" ) )
    pickle.dump(t,        open( "data_isotropic/tfull"    + suffix, "wb" ) )