
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
import pickle
import time
import os.path

########################################################################
# prints duration of time step and estimates completion time

def print_progress(t, tindex, t0, t1):
    time_of_step   = t1-t0
    num_steps_left = t.shape[0] - tindex
    est_time  = (time_of_step*num_steps_left/60.)/60.
    print('step took {0} seconds, about {1} hours remaining'.format(time_of_step, est_time)) 
    
########################################################################
# Generate random basepoints and timeline for initializing code
    
def blank_data(npoints, nparticles, t, x0, 
               savewhich = 'no history', which_database = 'channel'):
               
    if savewhich == 'history':
        x    = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles), dtype = np.float32)
        x[0] = x0
        LB   = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles), dtype = np.float32)
        LT   = np.zeros(shape = (t.shape[0] + 1, npoints, nparticles), dtype = np.float32)
    else:
        x  = x0
        LB = np.zeros(shape = (npoints, nparticles), dtype = np.float32)
        LT = np.zeros(shape = (npoints, nparticles), dtype = np.float32)
        
    disp   = np.zeros(shape = (t.shape[0], npoints   ), dtype = np.float32)
    HT     = np.ones( shape = (npoints, nparticles   ), dtype = np.float32)
    HT     = (2*t[0])*HT
    
    if which_database == 'channel':
        return x, HT, LB, LT, disp
    elif which_database == 'isotropic':
        return x, disp
        
########################################################################
# Save and Load data functions
        
def save_data(npoints, nparticles, kappa, tindex, t, x, LT, LB, HT, disp, savewhich = 'no history'):
              
    mainfolder = 'data_pure_brownian/'
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
          
    if not os.path.exists(folder):
        os.makedirs(folder)
        
    suffix = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, kappa) 
    pickle.dump(tindex,  open( folder + "tindex"  + suffix, "wb" ) )
    pickle.dump(t,       open( folder + "t"       + suffix, "wb" ) )
    pickle.dump(x,       open( folder + "x"       + suffix, "wb" ) )
    pickle.dump(LT,      open( folder + "LT"      + suffix, "wb" ) )
    pickle.dump(LB,      open( folder + "LB"      + suffix, "wb" ) )
    pickle.dump(HT,      open( folder + "HT"      + suffix, "wb" ) )
    pickle.dump(disp,    open( folder + "disp"    + suffix, "wb" ) )

    
def load_data(npoints, nparticles, kappa, savewhich = 'no history'):

    mainfolder = 'data_pure_brownian/'
        
    if savewhich == 'history':
        folder = mainfolder +'history/'
    else:
        folder = mainfolder
        
    suffix     = '_Points_{0}_Traj_{1}_Pr_{2}.p'.format(npoints, nparticles, kappa)
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
         
########################################################################
# Primary evolution function

def evolution(npoints, nparticles, kappa, t, x0, savewhich = 'no history'):
             
    mainfolder = 'data_pure_brownian/'
        
    if savewhich == 'history':
        folder = mainfolder + 'history/'
    else:
        folder = mainfolder
        
    suffix = '_Points_{0}_Traj_{1}_k_{2}.p'.format(npoints, nparticles, kappa)
    fname  =  folder + "tindex" + suffix
    if os.path.isfile(fname) == True:
        tindex, t, x, LT, LB, HT, disp  = load_data(npoints, nparticles, kappa, savewhich)
                                                        
        if (tindex+1) == disp.shape[0]:
            print 'Run Already Completed'
            return
        else:
            print 'Resuming Particle Simulations'
            evolve(kappa, tindex, t, x, HT, LB, LT, disp, savewhich)
    else:
        print 'Running Evolve from begining'
        x, HT, LB, LT, disp = blank_data(npoints, nparticles, t, x0, savewhich)
        evolve(kappa, 0, t, x, HT, LB, LT, disp, savewhich)       
        
                                      
########################################################################
# Do work for evolution and calculate physical quantities of interest

def Lepingle_alg(t, dt, x, dX, HT, LB, LT, kappa):
    '''This algorithm advances stochastic particles which stay inside a domain.
       It calculates the hitting times and Local times and appropriately modifies
       the trajectories so they can never pass top or bottom of the channel.     '''
       
    Bottom = 0 
    eps    = 0.0001

    # Record Hitting Time    
    cond1  = x[...] < Bottom + eps

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
    x[...] = np.clip(x[...], Bottom, 1000) 
    return x, HT, LB, LT
             
def Euler_Maruyama_step(t, dt, x, kappa = 1):
    dX  = 0
    dW  = np.random.randn(*x.shape)*np.sqrt(dt)
    dX += np.sqrt(2*kappa)*dW
    return x + dX, dX
             
def evolve_one_step(t, dt, x, HT, LB, LT, kappa = 1):

    x, dX         = Euler_Maruyama_step(t, 
                                        dt, 
                                        x, 
                                        kappa)
    x, HT, LB, LT = Lepingle_alg(t, 
                                 dt, 
                                 x, 
                                 dX, 
                                 HT, 
                                 LB, 
                                 LT, 
                                 kappa)
    return x, HT, LB, LT
    


def evolve(kappa, tind, t, x, HT, LB, LT, disp, savewhich = 'no history'): 
    dt    = abs(t[0] - t[1])
    if savewhich == 'history':
        npoints    = x.shape[1]
        nparticles = x.shape[2]
    else:
        npoints    = x.shape[0]
        nparticles = x.shape[1]
    
    T0 = time.time()  
    for tindex in range(tind, t.shape[0]):
        t0 = time.time()        
        print('Brownian step {0} of {1} for npoints {2}, nparticles {3} annd kappa Number = {4}'.format(tindex, t.shape[0], npoints, nparticles, kappa))   
        if savewhich == 'history':
            x[tindex + 1], HT, LB[tindex + 1], LT[tindex + 1] = evolve_one_step(t[tindex], 
                                                                                dt, 
                                                                                x[tindex], 
                                                                                HT, 
                                                                                LB[tindex], 
                                                                                LT[tindex], 
                                                                                kappa)                               
        else: 
            x, HT, LB, LT = evolve_one_step(t[tindex], 
                                            dt, 
                                            x, 
                                            HT, 
                                            LB, 
                                            LT, 
                                            kappa)            
        if tindex%100000 == 0:
            save_data(npoints, nparticles, kappa, tindex, t, x, LT, LB, HT, disp, savewhich)       
        t1 = time.time()
        print_progress(t, tindex, t0, t1) 
    T1 = time.time()
    print('FINISHED: kappa {0} took {1} hours total'.format(kappa, ((T1-T0)/60.)/60.)) 
               
    save_data(npoints, nparticles, kappa, tindex, t, x, LT, LB, HT, disp, savewhich) 