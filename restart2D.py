# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:05:57 2015

@author: glbjch
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import matplotlib
import os
import random
import pytoughgrav as ptg
import matplotlib.pyplot as plt
import time
import shutil

t0=time.clock()
plt.close('all')
#%%
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Steffi_GRAV')

mod='20150327_1_var'
basemod='20150327_1'
if not os.path.exists(mod):
    os.makedirs(mod)

# read template file    
dat=t2data(basemod+'/flow2.inp')
geo=mulgrid(basemod+'/grd.dat')
grid=dat.grid
ptg.makeradial(geo,None,width=10.)
# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
# Read from existing file
inc=t2incon(basemod+'/flow2.sav')
inc.write(mod+'/flow2.inc')

#%% update T2 params
dat.parameter['max_timestep']=2.0e6
dat.parameter['print_interval']=50
dat.parameter['timestep']=[1.0]
#dat.output_times['time_increment']
#dat.output_times['time']=[1.0]
#dat.output_times['num_times_specified']=len(dat.output_times['time'])
#dat.output_times['num_times']=7
dat.parameter['option'][12]=2

dat.parameter['print_block']='ay 40'#############################################################

#%% load recharge data
rech=np.loadtxt(r'C:\Users\glbjch\Local Documents\Work\Modelling\Steffi_GRAV\dev_files\norm_monthrech.txt')
## Define GENER block
#fpms=1.01308803322 # flux per meter squared
#%% evevation dependent params
fm=3.64742287695e-08
fc=2.15803271989e-06
#mingen=2.0e-7
#cols=[col for col in geo.columnlist]
#count=0
#%% run function
ptg.gen_variable(mod,geo,grid,dat,ts=rech,elev_m=fm,elev_c=fc)
#
#%% write files
geo.write(mod+'/grd.dat') 
#       
## write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk')
#   
## write tough2 input file   
dat.write(mod+'/flow2.inp')
shutil.copy('dev_files/initial_it2file',mod+'/'+mod)
print time.clock()-t0