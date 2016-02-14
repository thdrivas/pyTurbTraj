import numpy as np
import sympy as sp
import pickle

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import isotropic1024coarse as info

npoints = 2
nparticles = 1000
nsteps = info['time'].shape[0]

x = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
x[..., 0] = info['lx']*np.random.random(size = (npoints,))[:, None]
x[..., 1] = info['ynodes'][info['ynodes'].shape[0]//2]
x[..., 2] = info['lz']*np.random.random(size = (npoints,))[:, None]

kappa = (2*info['nu'])**.5

lJHTDB = libJHTDB()
lJHTDB.initialize()
subdivisions = 2
timeline = np.linspace(info['time'][-1], 0, num=subdivisions*nsteps+1)
xfull = np.zeros(shape = (subdivisions*nsteps+1, npoints, nparticles, 3),
                 dtype = np.float32)
xfull[0] = x
truncated = False
for tindex in range(subdivisions*nsteps):
    print('step {0}'.format(tindex))
    # get velocity
    u = lJHTDB.getData(
                timeline[tindex],
                xfull[tindex],
                sinterp = interpolation_code['M2Q8'],
                tinterp = interpolation_code['NoTInt'],
                data_set = info['name'],
                getFunction = 'getVelocity')
    # Euler Maruyama
    dt = timeline[tindex]-timeline[tindex+1]
    dW = np.random.randn(*xfull.shape[1:])*np.sqrt(dt)
    xfull[tindex+1] = xfull[tindex] - u*dt + kappa*dW
lJHTDB.finalize()

pickle.dump(xfull, open( "xfull.p", "wb" ) )
pickle.dump(timeline, open( "timeline.p", "wb" ) )