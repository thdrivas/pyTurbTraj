import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import channel as info

import pickle


nu = 5e-5
Prandtl = 1
kappa = nu/Prandtl
kolmogorov_time = 1e-3
kolmogorov_length = 1


Lx = info['lx']
Ly = info['ly']
Lz = info['lz']

Top = Ly/2
Bottom = -Ly/2

eps = Ly/100

npoints = 10
nparticles = 20
nsteps = info['time'].shape[0]
nsteps = nsteps


x         = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
x[..., 0] = Lx*np.random.random(size = (npoints,))[:, None]
x[..., 1] = info['ynodes'][info['ynodes'].shape[0]//10]
x[..., 2] = Lz*np.random.random(size = (npoints,))[:, None]


lJHTDB = libJHTDB()
lJHTDB.initialize()

T = info['time'][-1];
subdivisions = 2
t = np.linspace(info['time'][-1], info['time'][0], num=subdivisions*nsteps+1)
tau = t[0]-t
xfull = np.zeros(shape = (subdivisions*nsteps+1, npoints, nparticles, 3),
                 dtype = np.float32)
localtimeL = np.zeros(shape = (subdivisions*nsteps+1, npoints, nparticles), dtype = np.float32)
localtimeR = np.zeros(shape = (subdivisions*nsteps+1, npoints, nparticles), dtype = np.float32)

xfull[0] = x
for tindex in range(subdivisions*nsteps):
    print('step {0}'.format(tindex))
    # get velocity
    u = lJHTDB.getData(
                t[tindex],
                xfull[tindex],
                sinterp = interpolation_code['M2Q8'],
                tinterp = interpolation_code['NoTInt'],
                data_set = info['name'],
                getFunction = 'getVelocity')
    # Euler Maruyama with Reflection
    dt = t[tindex]-t[tindex+1]
    dW = np.random.randn(*xfull.shape[1:])*np.sqrt(dt)
    dX =  - u*dt + np.sqrt(2*kappa)*dW
    xfull[tindex+1] = xfull[tindex] + dX
    localtimeL[tindex+1] = localtimeL[tindex] 
    localtimeR[tindex+1] = localtimeR[tindex] 
    for k in range(npoints):
        for j in range(nparticles):
            if xfull[tindex][k,j,1] < Bottom + eps:
                VV = -2*dt*np.log(np.random.rand());
                Y = (-dX[k,j,:]+np.sqrt(2*kappa*VV+dX[k,j,:]**2))/2;
                dL = max(0, Y[1] - xfull[tindex][k,j,1] + Bottom);
                xfull[tindex+1][k,j,1] = xfull[tindex][k,j,1] + dL 
                localtimeL[tindex+1][k,j] = localtimeL[tindex+1][k,j] + dL/kappa
            if xfull[tindex][k,j,1] > Top - eps:
                VV = -2*dt*np.log(np.random.rand());
                Y = (dX[k,j,:]+np.sqrt(2*kappa*VV+dX[k,j,:]**2))/2;
                dL = max(0, Y[1] + xfull[tindex][k,j,1] - Top);
                xfull[tindex+1][k,j,1] = xfull[tindex][k,j,1] - dL 
                localtimeR[tindex+1][k,j] = localtimeR[tindex+1][k,j] + dL/kappa
            xfull[tindex+1][k,j,1]=max(min(xfull[tindex+1][k,j,1],Top),Bottom); 
lJHTDB.finalize()


pickle.dump(xfull, open( "Chan_xfull.p", "wb" ) )
pickle.dump(t, open( "Chan_t.p", "wb" ) )
pickle.dump(LL, open( "Chan_localtimeL.p", "wb" ) )
pickle.dump(LR, open( "Chan_localtimeR.p", "wb" ) )