
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
from tools import *
# functions that do the real work are defined in tools

#######################################################################


def main():     
    ''' Using this code, one can study molecular (with additive brownian noise)
        trajectories is turbulent flows using the JHU turbulence databases.  In
        particular, the code is currently set up to use the isotropic and channel 
        flow databases. We have used an algorithm of Lepingle to reflect particles
        at walls to modify the Euler Maruyama if studying particles in channel flow.
        
        This code will generate and save data on particle trajectories.  To analyse 
        data already generated, use the ipython notebooks analyze_statistics and
        analyze_history.                                                           '''
        
    npoints    = 2**2        # number of basepoints for particles
    nparticles = 2**9      # number of particle realizations at each basepoint

    # choose database of interest
    databases      = np.array(['channel', 'isotropic'])  
    which_database = databases[0]
    
    # select whether to store time-history of trajectories and local times 
    possible_saves = np.array(['history', 'no history']) 
    savewhich      = possible_saves[0]
    
    # randomly generate basepoints
    x0 = random_initial_x(npoints, 
                          nparticles, 
                          which_database)
    x0[..., 1] = -np.ones(x0[..., 1].shape)
    print x0
    
    # generates timeline, subdivides database timestep (velocity interpolation will be used)
    t  = get_timeline(which_database, 
                      subdiv = 3) 
    
    # array of Prandtl numbers
    PrandtlNumbers = np.array([1e1, 1e0, 1e-1])   
    for i in range(PrandtlNumbers.shape[0]):
        Prandtl = np.float(PrandtlNumbers[i])
        if check_if_complete(npoints, nparticles, Prandtl, savewhich, which_database) == 1:
            print which_database + ' Prandtl {0} already complete, moving to next'.format(Prandtl)
            continue       
        evolution(npoints, nparticles, Prandtl, t, x0, savewhich, which_database)
        
    compress_data(npoints, nparticles, PrandtlNumbers, t, savewhich, which_database)
    
if __name__ == "__main__":
    main()
    