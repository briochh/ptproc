# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\glbjch\.spyder2\.temp.py
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import os

mod='20140827_1_py_it'
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/')

if not os.path.exists(mod):
    os.makedirs(mod)

## top surface
#surf = np.loadtxt(
#r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\2Ddev\2dprof.txt',
#delimiter='\t', skiprows=1) # load surface file
surf=np.array([[0,500],[2000,500]])
min5=-5*np.ones((surf.shape[0],1)) # adapt to min max of y (-5,+5)
five=5*np.ones((surf.shape[0],1))
surf=np.concatenate(((np.concatenate((
np.hsplit(surf,2)[0],min5,np.hsplit(surf,2)[1]),axis=1)),
(np.concatenate((np.hsplit(surf,2)[0],five,np.hsplit(surf,2)[1]),axis=1))),
axis=0)

maxx=2000
origin=[0,0,510]

dx=25
dy=[10]
dz=25
nx=maxx/dx
zcells=[10]*34+[10]*25

geo = mulgrid().rectangular([dx]*nx,dy,zcells, origin=origin, atmos_type =0, 
convention = 2 )  # creates geometry 20 cells that are 500 m width in x,
# 1 cell 1000 m width in y
# 20 cells 100 m high in z  
# to make use of more possible grid names add char=ascii_lowercase+ascii_uppercase

geo.atmosphere_volume= 1.e50 # change volume of atmos cell to 1e50

geo.fit_surface(surf, silent=True, layer_snap=2.0) # fit topograpghy surface

# define and add wells
#well1=well('well1',[[50,5,725],[50,5,0]])
#well2=well('well2',[[500,5,600],[500,5,0]])
#well3=well('well3',[[1000,5,600],[1000,5,0]])
#well4=well('well4',[[1750,5,300],[1750,5,0]])
#well5=well('well5',[[2200,5,200],[2200,5,0]])
#well6=well('well6',[[2800,5,100],[2800,5,0]])
#geo.add_well(well1)
#geo.add_well(well2)
#geo.add_well(well3)
#geo.add_well(well4)
#geo.add_well(well5)
#geo.add_well(well6)

# write geometry to output file 
geo.write(mod+'/2dgrd.dat') 

###### MAKE TOUGH GRID
grid = t2grid().fromgeo(geo)

# define relative permeability and cp paramters to use
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}
norp={'type':5, 'parameters':[]}
cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.0]}
nocp={'type':1, 'parameters':[0.0,0.0,1.0]}

# define rock types and add cp and rp params
lp=rocktype('lp   ', nad=3, permeability = [1.e-16]*2+[1e-16],
porosity=0.1, conductivity=2.51, specific_heat=920) 
lp.dry_conductivity=1.5 
lp.tortuosity=0.0
lp.relative_permeability=rp
lp.capillarity=cp
grid.add_rocktype(lp)

hp=rocktype('hp   ', nad=3, permeability = [5.e-13]*2+[5.e-13],
porosity=0.34)
hp.dry_conductivity=1.5 
hp.tortuosity=0.0
hp.relative_permeability=rp
hp.capillarity=cp
grid.add_rocktype(hp)

b=rocktype('nocp ', nad=3, permeability = [5.e-13]*2+[5.e-13],
porosity=0.34)
b.dry_conductivity=1.5 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
grid.add_rocktype(b)

at=rocktype('atmos', nad=3, density=1.225, permeability = [5.e-13]*2+[5.e-13],
porosity=1.0)
at.dry_conductivity=1.5 
at.tortuosity=0.0
at.relative_permeability=norp
at.capillarity=nocp
grid.add_rocktype(at)

# assign rock properties
# define low permeability region
lam=0.004####################################################################
k0=5.0e-13
for blk in grid.blocklist[1:]:
#    if blk.centre[2] <= 250 and blk.centre[0] <= 1400: 
#       blk.rocktype = grid.rocktype['lp   ']
#    else: 
        blk.rocktype = grid.rocktype['hp   ']
        # permeability modification
        col=geo.column[geo.column_name(str(blk))]
        lay=geo.column_surface_layer(col)    
        hmax=geo.block_surface(lay,col)
        pmx=blk.rocktype.permeability[0]*np.exp(-lam*(hmax-blk.centre[2]))
        grid.block[(str(blk))].pmx=pmx
# define rocktype of atmospher block    
for blk in grid.atmosphere_blocks[:]:
    blk.rocktype= grid.rocktype['atmos']
    grid.block[(str(blk))].pmx=blk.rocktype.permeability[0]

# select last column in block list and set as no cp and rp. Set to large volume
bcol=geo.columnlist[-1]
for lay in geo.layerlist:
    blk=geo.block_name(lay.name, bcol.name)
    if blk in geo.block_name_list:
       grid.block[(blk)].rocktype= grid.rocktype['nocp ']
       grid.block[(blk)].volume=1E50
   
# read template file    
dat=t2data('initialflow2.inp')

# add rocktype, element and connection data to dat class
dat.grid=grid

# INCON
dat.incon.clear
# Define incon block
initP=1.013e5
initSG=0.99
initT=25.0
cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
dat.incon[geo.block_name_list[0]]=cond
for blk in grid.blocklist[1:]:
    if grid.block[str(blk)].rocktype==nocp:
       initP=1.013e5
       initSG=0.99
       initT=25.0
       cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
    elif blk.centre[2] < 0.0:
       initP=1.013e5+(997.0479*9.81*abs(blk.centre[2]))
       initSG=0.0
       initT=25.0
       cond=[[0.0,0.0,0.0],[initP,initSG,initT]]
       dat.incon[str(blk)]=cond
    elif grid.block[str(blk)].rocktype==lp:
       initP=1.013e5
       initSG=0.0
       initT=25.0
       cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond
    else:
       initP=1.013e5
       initSG=0.0
       initT=25.0
       cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond 
       
       

# Define GENER block
fpms=7.7354e-6 # flux per meter squared
fm=3.24e-8
fc=-7.199e-7
mingen=2.0e-7
cols=[col for col in geo.columnlist]
count=0

# time dependant generation
mult=0.9
yrsec=3600*24*365.25
sixmonth=yrsec/2
times=[0.0,1000*yrsec]+np.arange((1000*yrsec)+sixmonth,(1020*yrsec),sixmonth).tolist()+[(1020*yrsec),1.0e15]
numt=len(times)
dat.clear_generators()
for col in cols:
    count=count+1
    lay=geo.column_surface_layer(col)
    blkname=geo.block_name(lay.name,col.name)
    gx=(grid.block[blkname].centre[2]*fm)+fc
    if gx < mingen: gx=mingen# for elevation dependant recharge!
#    lowgx=(grid.block[blkname].centre[2]*fm)+(fc+(mult*fc))
#    if lowgx < mingen-(mult*mingen): lowgx=mingen-(mult*mingen)
#    highgx=(grid.block[blkname].centre[2]*fm)+(fc-(mult*fc))
#    if highgx < mingen+(mult*mingen): lowgx=mingen+(mult*mingen)
#    gxc=[gx]+((numt-3)/2)*[lowgx,highgx]+[gx,gx]
    ex=numt*[1.0942e5]
#    gxa=np.multiply(col.area,gxc).tolist()
#    gen=t2generator(name=' q'+col.name,block=blkname,type='COM1',gx=None,ex=None,hg=None,fg=None, rate=gxa, enthalpy=ex, time=times,ltab=numt,itab=numt-1)
    gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gx*col.area, ex=1.0942e5)
    dat.add_generator(gen)
    



       
# write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk',wells=True)
   
# write tough2 input file   
dat.write(mod+'/flow2.inp')
