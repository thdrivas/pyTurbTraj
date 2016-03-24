
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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from channel_script import *
from tools import *

import subprocess


########################################################################

npoints    = 10
nparticles = 100

databases      = np.array(['channel', 'isotropic'])  
which_database = databases[1]

if which_database == 'channel':
    DB = info
    file_title = '_Channel'
    folder = 'data_channel/history/'
    savefolder = 'figs/channel'
    moviefolder = 'movies/channel'
elif which_database == 'isotropic':
    DB = info_iso
    file_title = '_Isotropic'
    folder = 'data_isotropic/history/'
    savefolder = 'figs/isotropic'
    moviefolder = 'movies/isotropic'
    
Lx = DB['lx']
Ly = DB['ly']
Lz = DB['lz']

Top    = Ly/2
Bottom = -Ly/2

suffix  = file_title + '_data_points_{0}_particles_{1}.p'.format(npoints, nparticles)
    
x       = pickle.load( open( folder + "x"     + suffix, "rb" ) )
disp    = pickle.load( open( folder + "disp"  + suffix, "rb" ) )

PrandtlNumbers = np.array([1e-1, 1e0, 1e-1]) # also have 1e-3
colors = np.array(['Indianred', 'Steelblue','green', 'yellow'])  


########################################################################
# Make movies of puffs of particles evolving

for k in range(x[0].shape[1]):                   # base points 
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection = '3d')
    for i in range(x[0].shape[0]):               # times
        ax.cla() 
        ax.clear()
        for j in range(len(PrandtlNumbers)):     # Pr numbers
            ax.plot(   x[j][i, k, :, 0],
                       x[j][i, k, :, 1],
                       x[j][i, k, :, 2],
                       label = '$Pr = {0}$'.format(PrandtlNumbers[j]),
                       color = colors[j],
                       marker = 'o',
                       markersize = 5,
                       linestyle = 'None')
                       
        if which_database == 'channel':
            ax.set_xlim3d(-10 , 20)
            ax.set_ylim3d(-Ly/2, Ly/2)
            ax.set_zlim3d(6.5, 3*np.pi)
        elif which_database == 'isotropic':
            ax.set_xlim3d(-Lx/2, Lx/2)
            ax.set_ylim3d(-Ly/2, Ly/2)
            ax.set_zlim3d(-Lz/2, Lz/2)
    
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        #ax.view_init(elev=70., azim=-47)
        
        suffix = 'channel_npoints_{0}_nparticles_{1}_space_{2}'.format(npoints, nparticles, k)
        fig.savefig(savefolder + suffix + '_frame_{0}.png'.format(i+10), format = 'png')
    subprocess.call(['ffmpeg', '-framerate', '250', '-start_number', '11', '-i', savefolder + suffix + '_frame_%02d.png',
                     '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-r', '30', moviefolder + suffix + '.mp4'])