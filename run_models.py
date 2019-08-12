# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:35:03 2018

Script to run transport_func() many times!

@author: armitage
"""

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
from flem import initialise, solve_flem
import peakutils # https://zenodo.org/badge/latestdoi/102883046
import os

# Directory name
repatoir = 'executable_article/data/'
name = 'SFDctc'
dirtest = './%s/%s' % (repatoir, name)
wavename = '%s/wavelengths_%s.txt' % (repatoir,name)
sedname = '%s/sedfulxes_%s.txt' % (repatoir,name)
directory = os.path.dirname(dirtest)
if not os.path.exists(directory):
  os.makedirs(directory)

# Domain dimensions
dem = 0
bounds = [0, 2e5, 0, 8e5]
resolutions = [64, 128, 256, 512]

# Physical parameters
physical_space = [1e+0,1e-4,1.5,1,1e-4] # [kappa, c, nexp, alpha, U]

# 0 = MFD node-to-node; 1 = MFD cell-to-cell; 2 = SD node-to-node; 3 = SD cell-to-cell
flow = 3

# Time
dt = 1e4
num_steps = 501

# Output
out_time = 100

# Plot stuff
plot = 0

# Calculate valley to valley wavelength
statistics = 1

# output variables
number = np.linspace(0, 9, 10)
QS = np.zeros((len(resolutions),len(number),num_steps))
TIME = np.zeros((len(resolutions),len(number),num_steps))
WAVE = np.zeros((len(resolutions),len(number),20))

fwav = open(wavename,'w')
fsed = open(sedname,'w')

i = 0
for res in resolutions:
  j = 0
  for num in number:
    model_space, u_n, mesh, V, bc = initialise(dem, bounds, res)
    [QS[i,j,:],TIME[i,j,:],WAVE[i,j,:]] = solve_flem(model_space, physical_space, flow, u_n, mesh, V, bc, dt, num_steps, out_time, plot, statistics, name)
    k = 0
    for iw in range(20) :
      fwav.write('%d\t%d\t%g\n' % (res,num,WAVE[i,j,k]))
      k += 1
    j += 1
  i += 1

j = 0
for num in number :
  k = 0
  for res in resolutions:
    for t in range(num_steps):
      fsed.write('%d\t%d\t%g\t%g\n' % (num,res,TIME[k,j,t]*1e-6*ly*ly/kappa,QS[k,j,t]/dt*kappa))
    k += 1
  j += 1

fsed.close()

