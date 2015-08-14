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
import argparse

t0=time.clock()
plt.close('all')

#%%
parser = argparse.ArgumentParser(description='Prepare restart model')
parser.add_argument('-b','--base', help='basemodel name',required=True)
parser.add_argument('-l','--location', help='location',required=False, default='.')
parser.add_argument('-n','--number',help='var number',required=False,default='')
parser.add_argument('-t','--topsurf_flag', help='use pseudo top surface for recharge',action='store_true')
parser.add_argument('-e','--pseudo_elev', help='pseusdo constant elevation',required=False,default=None)
#parser.add_argument('-heat','--heatsource',help='heat source region for ptb',  nargs='+', type=int, required=False)
#parser.add_argument('-fr','--fluidsource',help='mass source injection?',required=False, default=False)
#parser.add_argument('-fr','--fluidsource',help='mass source injection?',required=False, default=False)



args = parser.parse_args()

#%%
os.chdir(args.location)


basemod=args.base
number=args.number
mod=basemod+'_var'+number
pseudo_topsurf=args.topsurf_flag
if not os.path.exists(mod):
    os.makedirs(mod)

# read template file    
dat=t2data('flow2.inp')
geo=mulgrid('grd.dat')
grid=dat.grid
ptg.makeradial(geo,None,width=10.)
# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
# Read from existing file
inc=t2incon('flow2.sav')
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
rech=np.loadtxt('20150527_2_rand.dat')
## Define GENER block
#fpms=1.01308803322 # flux per meter squared
#%% evevation dependent params
fm=3.64742287695e-08
fc=2.15803271989e-06
#mingen=2.0e-7
#cols=[col for col in geo.columnlist]
#count=0
#%% run function
#allgens,xs,zs,Areas,times=ptg.gen_variable(mod,geo,grid,dat,elev_m=fm,elev_c=fc,season_bias=0.7,new_rand=0.5)
if pseudo_topsurf:
    topsurf=np.loadtxt('2Dprofile.txt',delimiter='\t',skiprows=1)
    x=topsurf[:,0]
    z=topsurf[:,1]
    s=interpolate.UnivariateSpline(x,z)
    xnew=np.sort([col.centre[0] for col in geo.columnlist])
    znew=s(np.sort(xnew))
    plt.figure()
    plt.plot(x,z,xnew,znew)
    topsurf=np.vstack((xnew,znew)).T
else:
    topsurf=None

allgens,xs,zs,Areas,times=ptg.gen_variable(mod,geo,grid,dat,
                                           ts=rech,elev_m=fm,elev_c=fc,
                                           season_bias=0.7,pseudo_elev=float(args.pseudo_elev),pseudo_topsurf=None)
#
#%% write files
geo.write(mod+'/grd.dat') 
#       
## write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk')
#   
## write tough2 input file   
dat.write(mod+'/flow2.inp')
shutil.copy(basemod,mod+'/'+mod)
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

