# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 14:24:00 2015

@author: glbjch
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import matplotlib
import os
import random
import pytoughgrav as ptg
import ice_pytough as ipt
import matplotlib.pyplot as plt
import time
import shutil

t0=time.clock()
plt.close('all')
#%%
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi')


basemod='Coto20150911_1'
mod=basemod+'_ptb1'
if not os.path.exists(mod):
    os.makedirs(mod)
    
dat=t2data(basemod+'/flow2.inp')
geo=mulgrid(basemod+'/grd.dat')
grid=dat.grid
width=geo.bounds[1][1]-geo.bounds[0][1]
#ptg.makeradial(geo,None,width) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

yrsec=365.25*3600*24    
# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
inc=t2incon(basemod + '/flow2.sav')
for blk in inc.blocklist:
    dat.incon[blk]=[None,inc[blk][0:]]
dat.parameter['option'][12]= 0
dat.parameter['timestep']=[1.0, 1.0E3,8.6400e+04]
dat.output_times['time']=[1.0]
dat.output_times['time_increment']= 10*yrsec

dat.output_times['num_times_specified']=len(dat.output_times['time'])
dat.output_times['num_times']=200
dat.parameter['max_timestep']=1*yrsec # maximum timstep length
dat.parameter['print_interval']=50 # print (output) frequency to flow.out
dat.parameter['tstop']=5E3*yrsec

main=grid.rocktype['main ']
main.permeability=main.permeability*10

#for blk in grid.blocklist:
#    if blk.rocktype.name == 'main ':
#        #blk.pmx
#        blk.pmx=blk.pmx*10
        #blk.pmx

dat.clear_generators()
heat_flux=0.24
for blk in grid.blocklist[0:]: blk.hotcell=False
ipt.heatgen(mod,geo,dat,grid,heat_flux,function={'type':'log','points':[[5.0,1.],[10000.,0.24]]},inject=[150,1.0e-3,1.67e6])
ptg.gen_constant(mod,geo,grid,dat,constant=1.5e-5,enthalpy='var')#enthalpy=8440.)


#for col in geo.columnlist[:]:
#        if col not in ecol:
#    lay=geo.layerlist[-1] # bottom layer
#    blkname=geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
#    x=grid.block[blkname].centre[0]
#    if x <= 500:
#        cond=inc[blkname]
#        print blkname
#        initT=350
#        cond.variable[2]=initT
#        grid.block[blkname].volume=grid.block[blkname].volume*1E50

geo.write(mod + '/grd.dat')
grid.write_vtk(geo,mod+'/'+mod+'_initial.vtk') 
dat.incon.clear()
dat.write(mod + '/flow2.inp')

inc.write(mod + "/flow2.inc")


shutil.copy('dev_files/initial_it2file',mod+'/'+mod)
shutil.copy('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+basemod+'/genertot.txt',mod+'/genertot.txt')

