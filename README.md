# pyTurbTraj
numerically computes stochastic lagrangian particle trajectories in JHU turbulence databases

The script particles_script.p generates the data for given parameters.  
All of the functions that do the work called by this script are in tools.p.

After setting desired parameters in particles_script.p, run kickoff.py.  This script runs
particles_script.p and checks every 15 seconds to see if it's still running.  If it crashed due
to internet or database issues, it automatically restarts it.

Once the data saved, the data can be analyzed using the ipython notebooks analyze_history.ipynb 
and analyze_statistics.ipynb.
