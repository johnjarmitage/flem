from __future__ import print_function
import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from flem import initialise, solve_flem
import os

# Directory name
name = 'bergantes'
dirtest = './%s/' % name
directory = os.path.dirname(dirtest)
if not os.path.exists(directory) :
  os.makedirs(directory)


# Domain dimensions
dem = 1
bounds = [-0.539189, 40.432329, -0.019358, 41.007764]
res = 128
model_space, u_n, mesh, V, bc = initialise(dem, bounds, res)

# Physical parameters
physical_space = [1e-1,1e-5,1.5,0.1,1e-4] # [kappa, c, nexp, alpha, U]

# 0 = MFD node-to-node; 1 = MFD cell-to-cell; 2 = SD node-to-node; 3 = SD cell-to-cell
flow = 0

# Time
dt = 2e2
num_steps = 100

# Output
out_time = 100

# Plot stuff
plot = 0

# Calculate valley to valley wavelength
statistics = 0

[sed_flux, time, wavelength] = solve_flem(model_space,physical_space,flow,u_n, mesh, V, bc,dt,num_steps,out_time,plot,statistics,name)
