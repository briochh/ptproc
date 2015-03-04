## -*- coding: utf-8 -*-
#"""
#Spyder Editor
#
# Test script for running 1D conduction model. 
#"""
#
#
#
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import os # mudule for operating system comands

## Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
os.chdir(r"C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\Molly") # change directory
mod='20150202_1' # define model name
if not os.path.exists(mod): 
    os.makedirs(mod)    

perm=5.0e-14 # define permeability
poro=1e-7  # define porosity
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]} # relative permeability functions and parameters - if single phase not necessary  
norp={'type':5, 'parameters':[]} # no rel perm option

cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]} # capillary pressure functions and parameters - if single phase not necessary
nocp={'type':1, 'parameters':[0.0,0.0,1.0]} # no cp option

dx=1 # size of cell in x direction
dy=1 # size of cell in y direction
#origin=[0,0,750]
origin=[0,0,2000] # position of first cell in space
#zcells=[10]*20+[50]*10+[100]*6+[50]*10+[10]*20
zcells=[20]*100 # ([dimension] * number) of cells in z - can be changed to refine mesh in places    


## Define grid geometry object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rectangular
## atmos_type = 2 = no atmosphere 
## convention = 2 = 2 char for layer 3 digits for column - doesnt matter too much for small models
geo = mulgrid().rectangular([dx],[dy],zcells, origin=origin, atmos_type =2, 
convention = 2 ) 

# can be written to output file in model (mod) directory 
geo.write(mod+'/grd.dat') 

# MAKE TOUGH GRID 
grid = t2grid().fromgeo(geo)

## Create TOUGH input file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
dat=t2data('initialflow2.inp') # read from template file 
dat.parameter['print_block']=' e  1' # define element to print in output - useful for loggin progress of TOUGH sim 
dat.multi['num_equations']=3 # 3 defines non isothermal simulation

## define rock types - this just generates rock types they are not necessarily assigned to elements - that happens later 
rtypes=[] # creates empty list

# define object name for rock type. e.g. hp is the python object 'main ' will be the name of the ROCK in TOUGH input file. THE SPACE IN THE NAME IS IMPORTANT - MUST BE 5 char!
hp=rocktype('main ', nad=3, permeability = [perm]*2+[perm],
porosity=poro) 
hp.wet_conductivity=3 #  M/(m K) from Hickey - cotapaxi
hp.tortuosity=0.0
hp.relative_permeability=rp # if single phase this has no effect
hp.capillarity=cp # if single phase this has no effect
hp.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[hp] # add object to list of rocktypes

# define object for boundary rocktype - not used in current example.....
b=rocktype('bound', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
b.wet_conductivity=3 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
b.specific_heat=1000.0
rtypes=rtypes+[b]

# define object for base rock type - these rocktypes may all be the same but it can be useful to divide them up for defining initial conditions for different regions - there may be other ways to do this.....
base=rocktype('base ', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
base.wet_conductivity=3 
base.tortuosity=0.0
base.relative_permeability=norp
base.capillarity=nocp
base.specific_heat=1000.0
rtypes=rtypes+[base]

# define object for top rock type
top=rocktype('top  ', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
top.wet_conductivity=3 
top.tortuosity=0.0
top.relative_permeability=norp
top.capillarity=nocp
top.specific_heat=1000.0
rtypes=rtypes+[top]

# loop over list of rocktypes and add to pytough gird object
for rtype in rtypes:
    grid.add_rocktype(rtype)  

# Assign rock types and initial conditions to elements 
dat.incon.clear() # first clear any initial conditions that may have been defined before, for example in the template input file   
# lam=0.004 # can be used for permeability scaling with depth.
# k0=5.0e-13 # can be used for permeability scaling with depth.


if np.size(geo.columnlist) > 1: # can be used to find lateral boundaries in a 2D model - NOTE: WILL NOT WORK FOR 3D
    ecol=[geo.columnlist[0],geo.columnlist[-1]] # create list of boundary columns (first and last)
else: # if the column list length is only 1 then there can be no lateral boundary.
    ecol=[] # set boundary columns to none
    
# loop over every element block in the grid
for blk in grid.blocklist[0:]: 
    lay=geo.layer[geo.layer_name(str(blk))] # layer containing current block
    col=geo.column[geo.column_name(str(blk))] # column containing current block
    tlay=geo.column_surface_layer(col)    
    hmax=geo.block_surface(tlay,col)
#    if lay == geo.layerlist[1]
    if lay == tlay: # find surface elements    
        blk.rocktype=grid.rocktype['top  '] # assign rocktype "top  "
        grid.block[str(blk)].volume=1E50 # define as infinite VOLUME - THIS EFFECTIVELY FIXES PRESSURE AND TEMPERATURE FOR DURATION OF SIMULATION
        initP=5e4+(997.0479*9.81*abs(hmax-blk.centre[2])) # initial presure condition - MAY NOT BE APPROPRIATE - WHAT IS THE PRESSURE UNDER THICK GLACIER AT 5000 m amsl??
        initSG=0.0 # initial gas saturation   
        initT=1.0 # initial temperature - TOUGH2 doesn't seem to like < 1.0 C
        cond=[None,[initP,0.000,initT]] # initial condition object compiled    
    elif lay == geo.layerlist[-1]: # Find bottom layer
        blk.rocktype=grid.rocktype['base ']
        #grid.block[(blk)].volume=1E50
        initP=5e4+(997.0479*9.81*abs(hmax-blk.centre[2])) # initial pressure - function of depth...
        initSG=0.0
        initT=5.0
        #initT=0.85+((3.0/100.0)*abs(hmax-blk.centre[2])) # Can set up a variable temperature initial condition.
        cond=[None,[initP,0.000,initT]]        
    elif col in ecol: # if block is in lateral boundary columns
        blk.rocktype=grid.rocktype['bound']
        grid.block[str(blk)].volume=1E50
        initP=5e4+(997.0479*9.81*abs(hmax-blk.centre[2]))
        initSG=0.0
        initT=5.0
        #initT=0.85+((3.0/100.0)*abs(hmax-blk.centre[2]))
        cond=[None,[initP,0.000,initT]]
    else: # otherwise it is in the main model body    
        blk.rocktype = grid.rocktype['main ']
        initP=5e4+(997.0479*9.81*abs(hmax-blk.centre[2]))
        initSG=0.0
        initT=5.0
        #initT=0.85+((3.0/100.0)*abs(hmax-blk.centre[2]))
        cond=[None,[initP,0.000,initT]]
    dat.incon[str(blk)]=cond # add inital conditions for current block.
# for defining depth dependent permeability - JUST PERM - MAY WISH TO EXPLORE DEPTH DEPENDENT POROSITY 
#   pmx=blk.rocktype.permeability[0]*np.exp(-lam*(hmax-blk.centre[2])) # calculating depth dependent permeability modifier
#   grid.block[(str(blk))].pmx=pmx # Add PERMEABILITY MODIFIER for current element

# Add data stored in grid object to pyTOUGH dat object - this is what gets turned into the TOUGH input file        
dat.grid=grid    

#~ Create sources and sinks GENER block in TOUGH input file.      
dat.generator.clear() # first clear existing

tflux=0.09 # J/m2/s example heat flux

#cols=[col for col in geo.columnlist]
for col in geo.columnlist:
    lay=geo.layerlist[-1] # bottom layer
    blkname=geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
    gen=t2generator(name=' q'+col.name,block=blkname,type='HEAT',gx=tflux*col.area) # creat a generater oject with the heat generation rate of tflux - muliplication by column area important. 
    dat.add_generator(gen) # add generater to TOUGH2 input

# additional output parameters 
dat.parameter['max_timestep']=3.0e8 # maximum timstep length
dat.parameter['print_interval']=3 # print (output) frequency to flow.out
dat.parameter['timestep']=[1000.0] # initial timestep?
dat.output_times['time']=[1000.0,3600.0,8.6400e+04,3.1558e+07,3.1558e+08,3.1558e+09,3.1558e+10] # predefined output times
##dat.output_times['num_times_specified']=7
##dat.output_times['num_times']=7
#
#
#       
## write vtk of input information - can be view with para view
grid.write_vtk(geo,mod+'/inparam.vtk',wells=True) 
#   
## write tough2 input file   
dat.write(mod+'/flow2.inp')      
