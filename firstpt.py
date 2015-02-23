# -*- coding: utf-8 -*-
"""First model generation of TOUGH2 model with pytough"""
"""
Created on Fri May 16 15:52:34 2014

@author: glbjch
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *

## top surface
surf = np.loadtxt(
r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\2Ddev\2dprof.txt',
delimiter='\t', skiprows=1) # load surface file
min5=-5*np.ones((surf.shape[0],1)) # adapt to min max of y (-5,+5)
five=5*np.ones((surf.shape[0],1))
surf=np.concatenate(((np.concatenate((
np.hsplit(surf,2)[0],min5,np.hsplit(surf,2)[1]),axis=1)),
(np.concatenate((np.hsplit(surf,2)[0],five,np.hsplit(surf,2)[1]),axis=1))),
axis=0)

maxx=3500
origin=[0,0,750]

dx=10
dy=[10]
dz=10
nx=maxx/dx
zcells=[10]*34+[2]*80+[10]*19+[2]*30+[10]*25

geo = mulgrid().rectangular([dx]*nx,dy,zcells, origin=origin, atmos_type =0, 
convention = 2 )  # creates geometry 20 cells that are 500 m width in x,
# 1 cell 1000 m width in y
# 20 cells 100 m high in z  
# to make use of more possible grid names add char=ascii_lowercase+ascii_uppercase

geo.atmosphere_volume= 1.e50

# fit the sruface to box model
# geo.refine_layers(factor=2)
#geo.add_layer()
#refineZmin=250
#refineZmax=410
#
#lay=[geo.layer_index[geo.layer_containing_elevation(refineZmin).name],geo.layer_index[geo.layer_containing_elevation(refineZmax).name]]
#lay.sort()
#print lay
#reflaynum=np.concatenate(((np.arange(lay[0],lay[1])),[lay[1]]))
#
#reflay=[]
#for i in reflaynum:
##    toref=geo.layer_containing_elevation(z)
#    print i
#    reflay.append(geo.layer_name_from_number(i))
#    
#print reflay
##geo.refine_layers(layers=reflay,factor=10)
#    
#refineZmin=-10
#refineZmax=60
#
#lay=[geo.layer_index[geo.layer_containing_elevation(refineZmin).name],geo.layer_index[geo.layer_containing_elevation(refineZmax).name]]
#lay.sort()
#print lay
#reflaynum=np.concatenate(((np.arange(lay[0],lay[1])),[lay[1]]))
#
#reflay=[]
#for i in reflaynum:
##    toref=geo.layer_containing_elevation(z)
#    print i
##    reflay.append(geo.layer_name_from_number(i))
#
#reflay.append(geo.layer_name_from_number(reflaynum[-1]))    
#reflay.append(geo.layer_name_from_number(reflaynum[-2]))     
#print reflay
#geo.refine_layers(layers=reflay,factor=10)

geo.fit_surface(surf, silent=True, layer_snap=2.0)
well1=well('well1',[[50,5,725],[50,5,0]])
well2=well('well2',[[500,5,600],[500,5,0]])
well3=well('well3',[[1000,5,600],[1000,5,0]])
well4=well('well4',[[1750,5,300],[1750,5,0]])
well5=well('well5',[[2200,5,200],[2200,5,0]])
well6=well('well6',[[2800,5,100],[2800,5,0]])


geo.add_well(well1)
geo.add_well(well2)
geo.add_well(well3)
geo.add_well(well4)
geo.add_well(well5)
geo.add_well(well6)

geo.write('2dgrd.dat') # write geometry to output file 

grid = t2grid().fromgeo(geo)

rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,0.0,1.0]}
norp={'type':5, 'parameters':[]}
cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,0.0,0.0,0.0]}
nocp={'type':5, 'parameters':[]}

lp=rocktype('lp   ', nad=3, permeability = [1.e-16]*2+[1e-16],
porosity=0.1, conductivity=2.51, specific_heat=920) 
lp.dry_conductivity=1.5 
lp.tortuosity=0.0
lp.relative_permeability=rp
lp.capillarity=cp
grid.add_rocktype(lp)

hp=rocktype('hp   ', nad=3, permeability = [1.e-12]*2+[1.e-12],
porosity=0.34)
hp.dry_conductivity=1.5 
hp.tortuosity=0.0
hp.relative_permeability=rp
hp.capillarity=cp
grid.add_rocktype(hp)

b=rocktype('nocp ', nad=3, permeability = [1.e-12]*2+[1.e-12],
porosity=0.34)
b.dry_conductivity=1.5 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
grid.add_rocktype(b)

at=rocktype('atmos', nad=3, density=1.225, permeability = [1.e-12]*2+[1.e-12],
porosity=1.0)
at.dry_conductivity=1.5 
at.tortuosity=0.0
at.relative_permeability=norp
at.capillarity=nocp
grid.add_rocktype(at)

for blk in grid.blocklist[1:]:
    if blk.centre[2] <= 250 and blk.centre[0] <= 1400: 
       blk.rocktype = grid.rocktype['lp   ']
    else: blk.rocktype = grid.rocktype['hp   ']
for blk in grid.atmosphere_blocks[:]:
    blk.rocktype= grid.rocktype['atmos']
bcol=geo.columnlist[-1]
for lay in geo.layerlist:
    blk=geo.block_name(lay.name, bcol.name)
    if blk in geo.block_name_list:
       grid.block[(blk)].rocktype= grid.rocktype['nocp ']
       grid.block[(blk)].volume=1E50
   

    
dat=t2data('initialflow2.inp')
dat.grid=grid


dat.incon.clear

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
       initP=997.0479*9.81*abs(blk.centre[2])
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
       initSG=10.8
       initT=25.0
       cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
       dat.incon[str(blk)]=cond 

fpms=7.7354e-6 # flux per meter squared
cols=[col for col in geo.columnlist]
count=0
dat.clear_generators()
for col in cols:
    count=count+1
    lay=geo.column_surface_layer(col)
    blkname=geo.block_name(lay.name,col.name)
    gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=fpms*col.area, ex=1.0942e5)
    dat.add_generator(gen)
    



       

grid.write_vtk(geo,'inparam.vtk',wells=True)
   
dat.write('flow2.inp')



