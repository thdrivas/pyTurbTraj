
#######################################################################
#                                                                     #
#  Copyright 2016 Theodore D Drivas                                   #
#                                                                     #
#  This file is part of pyTurbTraj.                                   #
#                                                                     #
#  pyTurbTraj is free software: you can redistribute it and/or modify #
#  it under the terms of the GNU General Public License as published  #
#  by the Free Software Foundation, either version 3 of the License,  #
#  or (at your option) any later version.                             #
#                                                                     #
#  pyTurbTraj is distributed in the hope that it will be useful,      #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of     #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      #
#  GNU General Public License for more details.                       #
#                                                                     #
#  You should have received a copy of the GNU General Public License  #
#  along with pyTurbTraj.  If not, see <http://www.gnu.org/licenses/> #
#                                                                     #
#######################################################################



########################################################################
# Import Necessary Modules

import numpy as np
import sympy as sp

import pyJHTDB
from pyJHTDB import libJHTDB
from pyJHTDB.dbinfo import interpolation_code
from pyJHTDB.dbinfo import channel as info
from pyJHTDB.dbinfo import isotropic1024coarse as info_iso

import pickle
import time
import os.path

########################################################################
# prints duration of time step and estimates completion time

def print_progress(t, tindex, t0, t1):
    time_of_step   = t1-t0
    num_steps_left = t.shape[0] - tindex
    est_time_left  = (time_of_step*num_steps_left/60.)/60.
    print('took {0} seconds, about {1} hours remaining'.format(time_of_step, est_time_left)) 
    
########################################################################
# Generate random basepoints and timeline for initializing code
    
def random_initial_x(npoints, nparticles, which_database = 'channel'):
    
    if which_database == 'channel':
        DB = info
        N  = 10 # puts you near y = -1 wall
    elif which_database == 'isotropic':
        DB = info_iso
        N  = 2  # puts you in middle of box
        
    Lx = DB['lx']
    Ly = DB['ly']
    Lz = DB['lz']
    x0         = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
    x0[..., 0] = Lx*np.random.random(size = (npoints,))[:, None]
    x0[..., 1] = DB['ynodes'][DB['ynodes'].shape[0]//N]
    x0[..., 2] = Lz*np.random.random(size = (npoints,))[:, None]
    return x0

def get_timeline(which_database = 'channel', subdiv = 2):
    
    if which_database == 'channel':
        T      = 25.9935
        DB_dt  = 0.0065
        nsteps = int(T/DB_dt)
    elif which_database == 'isotropic':
        T      = info_iso['time'][-1];
        nsteps = info_iso['time'].shape[0] 
        
    return np.linspace(T, 0, num = subdiv*nsteps+1)

def blank_data(npoints, nparticles, t, x0, 
               savewhich = 'no history', which_database = 'channel'):
               
    if savewhich == 'history':
        x    = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles, 3), dtype = np.float32)
        x[0] = x0
        LB   = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles, 3), dtype = np.float32)
        LT   = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles, 3), dtype = np.float32)
    else:
        x  = x0
        LB = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
        LT = np.zeros(shape = (npoints, nparticles, 3), dtype = np.float32)
        
    disp   = np.zeros(shape = (t.shape[0], npoints   ), dtype = np.float32)
    HT     = np.ones( shape = (npoints, nparticles   ), dtype = np.float32)
    HT     = (2*t[0])*HT
    
    if which_database == 'channel':
        return x, HT, LB, LT, disp
    elif which_database == 'isotropic':
        return x, disp
        
########################################################################
# Save and Load data functions
        
def save_data(npoints, nparticles, Prandtl, tindex, t, x, LT, LB, HT, disp, 
              savewhich = 'no history'):
              
    mainfolder = 'data_channel/'
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
        
    suffix = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, Prandtl) 
    pickle.dump(tindex,  open( folder + "tindex"  + suffix, "wb" ) )
    pickle.dump(t,       open( folder + "t"       + suffix, "wb" ) )
    pickle.dump(x,       open( folder + "x"       + suffix, "wb" ) )
    pickle.dump(LT,      open( folder + "LT"      + suffix, "wb" ) )
    pickle.dump(LB,      open( folder + "LB"      + suffix, "wb" ) )
    pickle.dump(HT,      open( folder + "HT"      + suffix, "wb" ) )
    pickle.dump(disp,    open( folder + "disp"    + suffix, "wb" ) )
    
def save_data_iso(npoints, nparticles, Prandtl, tindex, t, x, disp, 
                  savewhich = 'no history'):
              
    mainfolder = 'data_isotropic/'
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
        
    suffix = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, Prandtl) 
    pickle.dump(tindex,  open( folder + "tindex"  + suffix, "wb" ) )
    pickle.dump(t,       open( folder + "t"       + suffix, "wb" ) )
    pickle.dump(x,       open( folder + "x"       + suffix, "wb" ) )
    pickle.dump(disp,    open( folder + "disp"    + suffix, "wb" ) )
    
def load_data(npoints, nparticles, Prandtl, 
              savewhich = 'no history', which_database = 'channel'):
    
    if which_database == 'channel':
        mainfolder = 'data_channel/'
    elif which_database == 'isotropic':
        mainfolder = 'data_isotropic/'
        
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
        
    suffix     = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, Prandtl)
    if which_database == 'channel':
        tindex     = pickle.load( open( folder + "tindex"   + suffix, "rb" ) )
        t          = pickle.load( open( folder + "t"        + suffix, "rb" ) )
        x          = pickle.load( open( folder + "x"        + suffix, "rb" ) )
        LT         = pickle.load( open( folder + "LT"       + suffix, "rb" ) )
        LB         = pickle.load( open( folder + "LB"       + suffix, "rb" ) )
        HT         = pickle.load( open( folder + "HT"       + suffix, "rb" ) )
        disp       = pickle.load( open( folder + "disp"     + suffix, "rb" ) )  
        return tindex, t, x, LT, LB, HT, disp  
    elif which_database == 'isotropic':
        tindex     = pickle.load( open( folder + "tindex"   + suffix, "rb" ) )
        t          = pickle.load( open( folder + "t"        + suffix, "rb" ) )
        x          = pickle.load( open( folder + "x"        + suffix, "rb" ) )
        disp       = pickle.load( open( folder + "disp"     + suffix, "rb" ) )  
        return tindex, t, x, disp
       
def check_if_complete(npoints, nparticles, Prandtl, 
                      savewhich = 'no history', which_database = 'channel'):
                      
    if which_database == 'channel':
        mainfolder = 'data_channel/'
    elif which_database == 'isotropic':
        mainfolder = 'data_isotropic/'
        
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
        
    suffix = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, Prandtl)
    fname  =  folder + "tindex" + suffix
    Bool = 0
    if os.path.isfile(fname) == True:
        if which_database == 'channel':
            tindex, t, x, LT, LB, HT, disp  = load_data(npoints, nparticles, Prandtl, savewhich, which_database)
        if which_database == 'isotropic':
            tindex, t, x, disp              = load_data(npoints, nparticles, Prandtl, savewhich, which_database)
        if (tindex+1) == disp.shape[0]:
            Bool = 1
    return Bool     
    
def compress_data(npoints, nparticles, PrandtlNumbers, t, 
                  savewhich = 'no history', which_database = 'channel'):
                  
    if which_database == 'channel':
        mainfolder = 'data_channel/'
    elif which_database == 'isotropic':
        mainfolder = 'data_isotropic/'

    if savewhich == 'history':
        folder = mainfolder + 'history/'
        if which_database == 'channel':
            dispf = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0],   npoints               ), dtype = np.float32)
            HTf   = np.zeros(shape = (PrandtlNumbers.shape[0],               npoints, nparticles   ), dtype = np.float32)
            Xf    = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0]+1, npoints, nparticles, 3), dtype = np.float32)
            LBf   = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0]+1, npoints, nparticles, 3), dtype = np.float32)
            LTf   = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0]+1, npoints, nparticles, 3), dtype = np.float32)
        elif which_database == 'isotropic':
            dispf = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0],   npoints               ), dtype = np.float32)
            Xf    = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0]+1, npoints, nparticles, 3), dtype = np.float32)
    else:
        folder = mainfolder
        if which_database == 'channel':
            dispf = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0], npoints               ), dtype = np.float32)
            HTf   = np.zeros(shape = (PrandtlNumbers.shape[0],             npoints, nparticles   ), dtype = np.float32)
            Xf    = np.zeros(shape = (PrandtlNumbers.shape[0],             npoints, nparticles, 3), dtype = np.float32)
            LBf   = np.zeros(shape = (PrandtlNumbers.shape[0],             npoints, nparticles, 3), dtype = np.float32)
            LTf   = np.zeros(shape = (PrandtlNumbers.shape[0],             npoints, nparticles, 3), dtype = np.float32)
        elif which_database == 'isotropic':
            dispf = np.zeros(shape = (PrandtlNumbers.shape[0], t.shape[0], npoints               ), dtype = np.float32)
            Xf    = np.zeros(shape = (PrandtlNumbers.shape[0],             npoints, nparticles, 3), dtype = np.float32)

    for i in range(PrandtlNumbers.shape[0]):
        if which_database == 'channel':        
            tindex, t, Xf[i], LTf[i], LBf[i], HTf[i], dispf[i]  = load_data(npoints, 
                                                                            nparticles, 
                                                                            PrandtlNumbers[i],
                                                                            savewhich, 
                                                                            which_database)
            suffix = '_Channel_data_points_{0}_particles_{1}.p'.format(npoints, nparticles)
            pickle.dump(dispf,  open( folder + "disp"  + suffix, "wb" ) )
            pickle.dump(Xf,     open( folder + "x"     + suffix, "wb" ) )
            pickle.dump(HTf,    open( folder + "HT"    + suffix, "wb" ) )
            pickle.dump(LBf,    open( folder + "LB"    + suffix, "wb" ) )
            pickle.dump(LTf,    open( folder + "LT"    + suffix, "wb" ) )                                                            
        elif which_database == 'isotropic':
            tindex, t, Xf[i], dispf[i]  = load_data(npoints, 
                                                    nparticles, 
                                                    PrandtlNumbers[i],
                                                    savewhich, 
                                                    which_database)
            suffix = '_Isotropic_data_points_{0}_particles_{1}.p'.format(npoints, nparticles)
            pickle.dump(dispf,  open( folder + "disp"  + suffix, "wb" ) )
            pickle.dump(Xf,     open( folder + "x"     + suffix, "wb" ) )              
                     
########################################################################
# Primary evolution function

def evolution(npoints, nparticles, Prandtl, t, x0, 
              savewhich = 'no history', which_database = 'channel'):
              
    if which_database == 'channel':
        mainfolder = 'data_channel/'
    elif which_database == 'isotropic':
        mainfolder = 'data_isotropic/'
        
    if savewhich == 'history':
        folder = mainfolder + 'history/'
    else:
        folder = mainfolder
        
    suffix = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, Prandtl)
    fname  =  folder + "tindex" + suffix
    if os.path.isfile(fname) == True:
        if which_database == 'channel':
            tindex, t, x, LT, LB, HT, disp  = load_data(npoints, nparticles, Prandtl,
                                                        savewhich, which_database)
        elif which_database == 'isotropic':
            tindex, t, x, disp              = load_data(npoints, nparticles, Prandtl, 
                                                        savewhich, which_database)
                                                        
        if (tindex+1) == disp.shape[0]:
            print 'Run Already Completed'
            return
        else:
            print 'Resuming Particle Simulations'
            if which_database == 'channel':
                evolve(Prandtl, tindex, t, x, HT, LB, LT, disp, savewhich)
            elif which_database == 'isotropic':
                evolve_nowall(Prandtl, tindex, t, x, disp, savewhich)
    else:
        print 'Running Evolve from begining'
        x, HT, LB, LT, disp = blank_data(npoints, nparticles, t, x0, savewhich)
        if which_database == 'channel':
           evolve(Prandtl, 0, t, x, HT, LB, LT, disp, savewhich)      
        elif which_database == 'isotropic':
           evolve_nowall(Prandtl, 0, t, x, disp, savewhich) 
        else:
            print 'choose a valid database'     
        
                                      
########################################################################
# Do work for evolution and calculate physical quantities of interest

def get_dispersion(x, npoints, nparticles):
    '''Calculates Mean Squared Dispersion'''
    numcombs   = np.float(nparticles*(nparticles-1)/2.)
    r = np.zeros(shape = (npoints, nparticles, nparticles), dtype = np.float32) 
    for k in range(npoints):
        for i in range(nparticles):
            for j in range(nparticles):
                if i<j:
                    r[k,i,j] = np.sum(np.square(x[k, i, :] - x[k, j, :]))           
    return np.sum(np.sum(r, axis=2),axis=1)/float(numcombs)

def Lepingle_alg(t, dt, x, dX, HT, LB, LT, kappa):
    '''This algorithm advances stochastic particles which stay inside a domain.
       It calculates the hitting times and Local times and appropriately modifies
       the trajectories so they can never pass top or bottom of the channel.     '''
       
    # Geometry of Channel
    Ly     =  info['ly']
    Top    =  Ly/2. #top of channel
    Bottom = -Ly/2. #bottom of channel
    eps    =  Ly/100.

    # Record Hitting Time    
    cond1  = x[..., 1] < Bottom + eps
    cond2  = x[..., 1] > Top    - eps
    T      = 25.9935
    tau    = T - t
    HT[np.logical_or(cond1, cond2)] = np.minimum(HT[np.logical_or(cond1, cond2)], tau)
    
    # Use Lepingle's Algorithm for Local time  
        #Lepingle, D 1995 Euler scheme for reflected stochastic differential equations. 
        #Mathematics and Computers in Simulation 38 (1), 119-126.
    c1indices = np.where(cond1)   
    if c1indices[0].size > 0:
        V                    =  -2*dt*np.log(np.random.random(x.shape))[c1indices]
        Y                    = (-dX[c1indices] + np.sqrt(2*kappa*V + dX[c1indices]**2))/2.0
        dL                   =  Y - x[c1indices] + Bottom
        dL[np.where(dL < 0)] = 0
        x[c1indices]        += dL
        LB[c1indices]       += dL/kappa
    c2indices = np.where(cond2) 
    if c2indices[0].size > 0:
        V                    =  -2*dt*np.log(np.random.random(x.shape))[c2indices]
        Y                    = (-dX[c2indices] + np.sqrt(2*kappa*V + dX[c2indices]**2))/2.0
        dL                   =  Y + x[c2indices] - Top
        dL[np.where(dL < 0)] = 0
        x[c2indices]        -= dL
        LT[c2indices]       += dL/kappa
    x[..., 1] = np.clip(x[..., 1], Bottom, Top) 
    return x, HT, LB, LT
             
def Euler_Maruyama_step(lJHTDB, t, dt, x, kappa = 1, 
                        evolve_backward = True, which_database = 'channel'):
    if which_database == 'channel':
        DB = info
    elif which_database == 'isotropic':
        DB = info_iso
    '''Takes an Euler-Maruyama time step, backwards or forwards with additive noise'''
    u = lJHTDB.getData(
                t,
                x,
                sinterp     = interpolation_code['M2Q8'],
                tinterp     = interpolation_code['NoTInt'],
                data_set    = DB['name'],
                getFunction = 'getVelocity')
    if evolve_backward:
        dX = -u*dt
    else:
        dX = u*dt
    dW  = np.random.randn(*x.shape)*np.sqrt(dt)
    dX += np.sqrt(2*kappa)*dW
    return x + dX, dX
             
def evolve_one_step(lJHTDB, t, dt, x, HT, LB, LT, kappa = 1, 
                    evolve_backward = True):

    x, dX         = Euler_Maruyama_step(lJHTDB,
                                        t, 
                                        dt, 
                                        x, 
                                        kappa, 
                                        evolve_backward,
                                        'channel')
    x, HT, LB, LT = Lepingle_alg(t, 
                                 dt, 
                                 x, 
                                 dX, 
                                 HT, 
                                 LB, 
                                 LT, 
                                 kappa)
    return x, HT, LB, LT
    
    
def evolve_one_step_iso(lJHTDB, t, dt, x, kappa = 1, 
                        evolve_backward = True):

    x, dX         = Euler_Maruyama_step(lJHTDB,
                                        t, 
                                        dt, 
                                        x, 
                                        kappa, 
                                        evolve_backward,
                                        'isotropic')
    return x



def evolve(Prandtl, tind, t, x, HT, LB, LT, disp, 
           savewhich = 'no history', evolve_backward = True): 
    nu    = 5e-5   
    kappa = nu/Prandtl
    dt    = abs(t[0] - t[1])
    if savewhich == 'history':
        npoints    = x.shape[1]
        nparticles = x.shape[2]
    else:
        npoints    = x.shape[0]
        nparticles = x.shape[1]
    
    lJHTDB = libJHTDB()
    lJHTDB.initialize()
    
    T0 = time.time()  
    for tindex in range(tind, t.shape[0]):
        t0 = time.time()        
        print('step {0} of {1} for npoints {2}, nparticles {3} annd Prandtl Number = {4}'.format(tindex, t.shape[0], npoints, nparticles, Prandtl))   
        if savewhich == 'history':
            x[tindex + 1], HT, LB[tindex + 1], LT[tindex + 1] = evolve_one_step(lJHTDB,
                                                                                t[tindex], 
                                                                                dt, 
                                                                                x[tindex], 
                                                                                HT, 
                                                                                LB[tindex], 
                                                                                LT[tindex], 
                                                                                kappa, 
                                                                                evolve_backward)                               
            disp[tindex]  = get_dispersion(x[tindex], npoints, nparticles)
        else: 
            x, HT, LB, LT = evolve_one_step(lJHTDB,
                                            t[tindex], 
                                            dt, 
                                            x, 
                                            HT, 
                                            LB, 
                                            LT, 
                                            kappa, 
                                            evolve_backward)            
            disp[tindex]  = get_dispersion(x, npoints, nparticles)
        if tindex%100 == 0:
            save_data(npoints, nparticles, Prandtl, tindex, t, x, LT, LB, HT, disp, savewhich)       
        t1 = time.time()
        print_progress(t, tindex, t0, t1) 
    T1 = time.time()
    print('FINISHED: Prandtl {0} took {1} hours total'.format(Prandtl, ((T1-T0)/60.)/60.)) 
               
    lJHTDB.finalize() 
    save_data(npoints, nparticles, Prandtl, tindex, t, x, LT, LB, HT, disp, savewhich) 
    
def evolve_nowall(Prandtl, tind, t, x, disp, 
               savewhich = 'no history', evolve_backward = True): 
    nu    = info_iso['nu']
    kappa = nu/Prandtl
    dt    = abs(t[0] - t[1])
    if savewhich == 'history':
        npoints    = x.shape[1]
        nparticles = x.shape[2]
    else:
        npoints    = x.shape[0]
        nparticles = x.shape[1]
    
    lJHTDB = libJHTDB()
    lJHTDB.initialize()
    
    T0 = time.time()  
    for tindex in range(tind, t.shape[0]):
        t0 = time.time()        
        print('step {0} of {1} for npoints {2}, nparticles {3} annd Prandtl Number = {4}'.format(tindex, t.shape[0], npoints, nparticles, Prandtl))      
        if savewhich == 'history':
            x[tindex + 1] = evolve_one_step_iso(lJHTDB,
                                                t[tindex], 
                                                dt, 
                                                x[tindex], 
                                                kappa, 
                                                evolve_backward)                               
            disp[tindex]  = get_dispersion(x[tindex], npoints, nparticles)
        else: 
            x = evolve_one_step_iso(lJHTDB,
                                    t[tindex], 
                                    dt, 
                                    x, 
                                    kappa, 
                                    evolve_backward)            
            disp[tindex]  = get_dispersion(x, npoints, nparticles)
        if tindex%100 == 0:
            save_data_iso(npoints, nparticles, Prandtl, tindex, t, x, disp, savewhich)       
        t1 = time.time()
        print_progress(t, tindex, t0, t1) 
    T1 = time.time()
    print('FINISHED: Prandtl {0} took {1} hours total'.format(Prandtl, ((T1-T0)/60.)/60.)) 
               
    lJHTDB.finalize() 
    save_data_iso(npoints, nparticles, Prandtl, tindex, t, x, disp, savewhich)   
      