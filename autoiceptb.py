# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:43:04 2015

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
import argparse

t0=time.clock()
plt.close('all')
#%%

parser = argparse.ArgumentParser(description='Prepare perturbation model')
parser.add_argument('-b','--base', help='basemodel name',required=True)
parser.add_argument('-l','--location', help='location',required=False, default='.')

args = parser.parse_args()

os.chdir(args.location)


basemod=args.base
mod=basemod+'_ptb'
if not os.path.exists(mod):
    os.makedirs(mod)
    
dat=t2data(basemod+'/flow2.inp')
geo=mulgrid(basemod+'/grd.dat')
grid=dat.grid
width=geo.bounds[1][1]-geo.bounds[0][1]
ptg.makeradial(geo,None,width)

yrsec=365.25*3600*24    
# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
inc=t2incon(basemod + '/flow2.sav')

dat.parameter['option'][12]= 0
dat.parameter['timestep']=[1.0, 1.0E3,8.6400e+04]
dat.output_times['time']=[1.0]
dat.output_times['time_increment']= 10*yrsec

dat.output_times['num_times_specified']=1
dat.output_times['num_times']=200
dat.parameter['max_timestep']=10*yrsec # maximum timstep length
dat.parameter['print_interval']=50 # print (output) frequency to flow.out
dat.parameter['tstop']=1E3*yrsec


for col in geo.columnlist[:]:
#        if col not in ecol:
    lay=geo.layerlist[-1] # bottom layer
    blkname=geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
    x=grid.block[blkname].centre[0]
    if x <= 500:
        cond=inc[blkname]
        print blkname
        initT=350
        cond.variable[2]=initT
        grid.block[blkname].volume=grid.block[blkname].volume*1E50

geo.write(mod + '/grd.dat')
grid.write_vtk(geo,mod+'/'+mod+'_initial.vtk') 
dat.write(mod + '/flow2.inp')

inc.write(mod + "/flow2.inc")


shutil.copy(basemod+'/'+basemod,mod+'/'+mod)
shutil.copy(basemod+'/genertot.txt',mod+'/genertot.txt')

