# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 15:41:29 2015

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
#os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi')

parser = argparse.ArgumentParser(description='Prepare perturbation model')
parser.add_argument('-b','--base', help='basemodel name',required=True)
parser.add_argument('-l','--location', help='location',required=False, default='.')
parser.add_argument('-heat','--heatsource',help='heat source region for ptb',required=False, default=[[0,0,3000],[500,0,3050]])

args = parser.parse_args()
heatsource=args.heatsource
os.chdir(args.location)

basemod=args.base+'_ptb'
#basemod='Cota20150612_1_ptb'
mod=args.base+'_rtn'
if not os.path.exists(mod):
    os.makedirs(mod)
    
dat=t2data('flow2.inp')
geo=mulgrid('grd.dat')
grid=dat.grid
width=geo.bounds[1][1]-geo.bounds[0][1]
ptg.makeradial(geo,None,width)

yrsec=365.25*3600*24    
# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
inc=t2incon('flow2.sav')

# additional output parameters 
dat.parameter['max_timestep']=3.1558e+09 # maximum timstep length
dat.parameter['print_interval']=20 # print (output) frequency to flow.out
dat.parameter['timestep']=[1.0E3,8.6400e+04,3.1558e+06] # initial timestep?
dat.parameter['tstop']=None
dat.output_times['time']=[1000.0,3600.0,8.6400e+04,3.1558e+07] # predefined output times
dat.output_times['num_times_specified']=len(dat.output_times['time'])
dat.output_times['num_times']=150
dat.output_times['time_increment']= 100*yrsec

for blk in grid.blocklist[0:]:
    lay=geo.layer[geo.layer_name(str(blk))] # layer containing current block
    col=geo.column[geo.column_name(str(blk))] # column containing current block
    blkname=blk.name #geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
    if (heatsource is not None and 
        blk.centre[2] > heatsource[0][2] and 
        blk.centre[2] <= heatsource[1][2] and 
        blk.centre[0] > heatsource[0][0] and 
        blk.centre[0] <= heatsource[1][0]): # if in heatsource region
        print blkname
        blk.volume=blk.volume/1E50

geo.write(mod + '/grd.dat')
grid.write_vtk(geo,mod+'/'+mod+'_initial.vtk') 
dat.write(mod + '/flow2.inp')

inc.write(mod + "/flow2.inc")

shutil.copy(basemod,mod+'/'+mod)
shutil.copy('genertot.txt',mod+'/')
