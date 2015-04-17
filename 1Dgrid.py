# -*- coding: utf-8 -*-
"""Function to create 1D TOUGH2 test model with pytough"""
"""
Created on Fri Jun 06 14:20:48 2014

@author: glbjch
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import os
import shutil

os.chdir('/Users/briochh/Documents/Workhere/testing')
mod='1Dtest3'
if not os.path.exists(mod):
    os.makedirs(mod)

perm=5.0e-13
poro=0.34
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}
cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]}
rechshift=0
wt=50.0
dx=1
dy=1
#origin=[0,0,750]
origin=[0,0,145]
zcells=[5]+[1]*140 
#zcells=[10]*34+[2]*80+[10]*19+[2]*30+[10]*25    

geo = mulgrid().rectangular([dx],[dy],zcells, origin=origin, atmos_type =0, 
convention = 2)
surface=np.array([[   0. ,   -0.5,  143. ],
       [   1. ,   -0.5,  143. ],
       [   0. ,    0.5,  143. ],
       [   1. ,    0.5,  143. ]])
geo.fit_surface(surface, silent=True, layer_snap=2.0) # fit topograpghy surface

geo.atmosphere_volume= 1.0e50

# write geometry to output file 
geo.write(mod+'/grd.dat') 

###### MAKE TOUGH GRID
grid = t2grid().fromgeo(geo)

# define relative permeability and cp paramters to use
#rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}
norp={'type':5, 'parameters':[]}
#cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]}
nocp={'type':1, 'parameters':[0.0,0.0,1.0]}



# define rock types and add cp and rp params
lp=rocktype('lp   ', nad=3, permeability = [1.e-16]*2+[1e-16],
porosity=0.1, conductivity=2.51, specific_heat=920) 
lp.dry_conductivity=1.5 
lp.tortuosity=0.0
lp.relative_permeability=rp
lp.capillarity=cp
grid.add_rocktype(lp)

hp=rocktype('hp   ', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
hp.dry_conductivity=1.5 
hp.tortuosity=0.0
hp.relative_permeability=rp
hp.capillarity=cp
grid.add_rocktype(hp)

b=rocktype('nocp ', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
b.dry_conductivity=1.5 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
grid.add_rocktype(b)

at=rocktype('atmos', nad=3, density=1.225, permeability = [perm]*2+[perm],
porosity=1.0)
at.dry_conductivity=1.5 
at.tortuosity=0.0
at.relative_permeability=norp
at.capillarity=nocp
grid.add_rocktype(at)


# define rocktype of atmospher block    
for blk in grid.atmosphere_blocks[:]:
    blk.rocktype= grid.rocktype['atmos']
    grid.block[(str(blk))].pmx=blk.rocktype.permeability[0]
    
    # assign rock properties
# define low permeability region
lam=0.004
k0=5.0e-13
for blk in grid.blocklist[1:]:
    blk.rocktype = grid.rocktype['hp   ']
        # permeability modification
    col=geo.column[geo.column_name(str(blk))]
    lay=geo.column_surface_layer(col)    
    hmax=geo.block_surface(lay,col)
    pmx=blk.rocktype.permeability[0]*np.exp(-lam*(hmax-blk.centre[2]))
    grid.block[(str(blk))].pmx=pmx
    
blay=geo.layerlist[-1]
for col in geo.columnlist:
    blk=geo.block_name(blay.name,col.name)
    if blk in geo.block_name_list:
       grid.block[(blk)].volume=1E50
    
# read template file    
dat=t2data('initialflow2.inp')
dat.parameter['print_block']='ee  1'
# add rocktype, element and connection data to dat class
dat.grid=grid    

# INCON
dat.incon.clear()
# Define incon block
initP=1.013e5
initSG=0.99
initT=25.0
cond=[None,[1.013e5,initSG,initT]]
dat.incon[geo.block_name_list[0]]=cond
for blk in grid.blocklist[1:]:
    if grid.block[str(blk)].rocktype==nocp:
       initP=1.013e5
       initSG=0.99
       initT=25.0
       cond=[None,[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
    elif blk.centre[2] < 50.0:
       initP=1.013e5+(997.0479*9.81*abs(50-blk.centre[2]))
       initSG=0.0
       initT=25.0
       cond=[None,[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
    elif grid.block[str(blk)].rocktype==lp:
       initP=1.013e5
       initSG=0.0
       initT=25.0
       cond=[None,[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
    else:
       initP=1.013e5
       initSG=0.0
       initT=25.0
       cond=[None,[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
       
dat.generator.clear()

# Define GENER block
fpms=7.7354e-6 # flux per meter squared
fm=3.24e-8
fc=-7.199e-7


mingen=2.0e-7
cols=[col for col in geo.columnlist]
count=0

dat.clear_generators()
for col in cols:
    count=count+1
    lay=geo.column_surface_layer(col)
    blkname=geo.block_name(lay.name,col.name)
    #gx=(grid.block[blkname].centre[2]*fm)+fc
    gx=fpms
    if gx < mingen: gx=mingen# for elevation dependant recharge!
    ex=1.0942e5
    gen=t2generator(name=' q'+col.name,block=blkname,type='COM1',gx=gx*col.area,ex=ex,hg=None,fg=None)
    #gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gx*col.area, ex=1.0942e5)
    dat.add_generator(gen) 

dat.parameter['max_timestep']=3.0e6
dat.parameter['print_interval']=30
#dat.parameter['timestep']=[1000.0]
#dat.output_times['time']=[1000.0,3600.0,8.6400e+04,3.1558e+07,3.1558e+08,3.1558e+09,3.1558e+10]
#dat.output_times['num_times_specified']=7
#dat.output_times['num_times']=7


       
# write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk',wells=True)
   
# write tough2 input file   
dat.write(mod+'/flow2.inp')      
shutil.copy('itough_input',mod+'/'+mod)