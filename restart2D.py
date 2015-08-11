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
from scipy import interpolate

t0=time.clock()
plt.close('all')
#%%
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper')

mod='20150806_4_var2'
basemod='20150806_4'
psuedo_topsurf=False
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
dat.output_times['time_increment']=3*28.*3600*24
dat.output_times['time']=[1.0]
dat.output_times['num_times_specified']=len(dat.output_times['time'])
dat.output_times['num_times']=500
dat.parameter['option'][12]=2
dat.parameter['tstop']=1.5E2*3600*24*365.32

dat.parameter['print_block']='ay 40'#############################################################

#%% load recharge data
rech=np.loadtxt(r'C:\Users\glbjch\Local Documents\Work\Modelling\Gravpaper\dev_files\20150527_2_rand.dat')
## Define GENER block
#fpms=1.01308803322 # flux per meter squared
#%% evevation dependent params
fm=3.64742287695e-08
fc=2.15803271989e-06
#mingen=2.0e-7
#cols=[col for col in geo.columnlist]
#count=0
#%% run function

if psuedo_topsurf:
    topsurf=np.loadtxt('dev_files/2Dprofile.txt',delimiter='\t',skiprows=1)
    x=topsurf[:,0]
    z=topsurf[:,1]
    s=interpolate.UnivariateSpline(x,z)
    xnew=np.sort([col.centre[0] for col in geo.columnlist])
    znew=s(np.sort(xnew))
    plt.figure()
    plt.plot(x,z,xnew,znew)
    topsurf=np.vstack((xnew,znew)).T

#allgens,xs,zs,Areas,times=ptg.gen_variable(mod,geo,grid,dat,elev_m=fm,elev_c=fc,season_bias=0.7,new_rand=0.5)
allgens,xs,zs,Areas,times=ptg.gen_variable(mod,geo,grid,dat,
                                           ts=rech,elev_m=fm,elev_c=fc,
                                           season_bias=0.7,psuedo_elev=400.0)
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

allgenscaled=np.divide(allgens.T,Areas)
inds=np.array(xs).argsort()
scaled_sorted=allgenscaled.T[inds]
fig,ax1=plt.subplots()
im=ax1.pcolormesh(np.sort(xs),np.array(times)/3600/24/365.25,scaled_sorted.T,vmin=np.min(scaled_sorted),cmap='rainbow')
cbar=fig.colorbar(im,ax=ax1,format="%.1e")
ax1.set_ylabel
ax1.set_xlabel
ax1.set_ylim((0,110))
ax1.set_xlim((0,np.max(xs)))
ax2=plt.twinx(ax1)
ax2.scatter(xs,zs)
ax2.set_xlim((0,np.max(xs)))

fig.savefig(mod+'/Variable_recharge.pdf')

