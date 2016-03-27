
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

import os
import time
import subprocess
import sys

#######################################################################                
        
def main():
    
    # This code starts and checks on the status of the main script, particles_script.py
    # Occationally, JHUTDB goes down momentarily either on the users end (internet lags)
    # or on the database side.  This script checks for the "heartbeat" of the code and,
    # if for whatever reason it died, will automatically restart it.    
    
	proc = restart_process()
	while True:
		time.sleep(15)

		if proc.poll() is None:
			print "Process is ok!"
		else:
			print "Process died! Restarting it." 
			proc = restart_process()


def restart_process():
	proc = subprocess.Popen([sys.executable, "./particles_script.py", r"&"], stdout=sys.stdout)
	return proc


if __name__ == "__main__":
	main()