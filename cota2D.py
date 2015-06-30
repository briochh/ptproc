





"""
Created on Mon Mar  9 14:21:33 2015

@author: molly
"""

import ice_pytough as ipt
import os
from t2grids import *
from t2data import *
import time 
import shutil
import pytoughgrav as ptg

#%% Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0=time.clock()
os.chdir("C:\Users\glbjch\Local Documents\Work\Modelling\Cotapaxi") # define working directory
mod='Cota20150630_1'
if not os.path.exists(mod):
    os.makedirs(mod)
#%%
yrsec=3600*365.25*24
origin=[0,0,6000] # position of first cell in space
width=1.0
zcells=[10]*35+[5]*5+[2]*20+[5]*5+[10]*86+[50]*4+[100]*15    #+[100]*6+[50]*10+[10]*20    
dy=1 # size of cell in y direction
xcells=[5]*12+[10]*254+[50]*27+[100,200,300,400,500,600,700,800,900,1000]  #+[100]*6+[50]*10+[10]*20

surf=ptg.topsurf('dev_files/Topography_crater_f.txt',delim='\t',headerlines=1,width=width)

geo=ipt.icegeo( mod, width=width, celldim=10., origin=origin,
              zcells=zcells, xcells=xcells, surface=surf, atmos_type=1, min_thick=2.0)

#%%
dat=t2data('dev_files/initialflow2.inp') # read from template file 
dat.parameter['print_block']=' w 46' # define element to print in output - useful for loggin progress of TOUGH sim 
dat.multi['num_equations']=3 # 3 defines non isothermal simulation

perm=1.0e-13 # define permeability
poro=0.1  # define porosity
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]} # relative permeability functions and parameters - if single phase not necessary  
norp={'type':5, 'parameters':[]} # no rel perm option
cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]} # capillary pressure functions and parameters - if single phase not necessary
nocp={'type':1, 'parameters':[0.0,0.0,1.0]} # no cp option

conds=4.0
heat_flux=0.5 #0.24

# MAKE TOUGH GRID 


#%% define rock types - this just generates rock types they are not necessarily assigned to elements - that happens later 
rtypes=[] # creates empty list
dat.incon.clear()
# define object name for rock type. e.g. hp is the python object 'main ' will be the name of the ROCK in TOUGH input file. THE SPACE IN THE NAME IS IMPORTANT - MUST BE 5 char!
main=rocktype('main ', nad=3, permeability = [perm]*2+[perm],
porosity=poro) 
main.conductivity= 4 #  M/(m K) from Hickey - cotapaxi
main.tortuosity=0.0
main.relative_permeability=rp # if single phase this has no effect
main.capillarity=cp # if single phase this has no effect
main.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[main] # add object to list of rocktypes

main=rocktype('surf ', nad=3, permeability = [perm]*2+[perm],
porosity=poro) 
main.conductivity= 4 #  M/(m K) from Hickey - cotapaxi
main.tortuosity=0.0
main.relative_permeability=rp # if single phase this has no effect
main.capillarity=cp # if single phase this has no effect
main.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[main] # add object to list of rocktypes

# define object for boundary rocktype - not used in current example.....
b=rocktype('bound', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
b.conductivity=4 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
b.specific_heat=1000.0
rtypes=rtypes+[b]

# define object for base rock type - these rocktypes may all be the same but it can be useful to divide them up for defining initial conditions for different regions - there may be other ways to do this.....
source=rocktype('sourc', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
source.conductivity=4 
source.tortuosity=0.0
source.relative_permeability=norp
source.capillarity=nocp
source.specific_heat=1000.0
rtypes=rtypes+[source]

hotcell=rocktype('hotcl', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
hotcell.conductivity=4 
hotcell.tortuosity=0.0
hotcell.relative_permeability=norp
hotcell.capillarity=nocp
hotcell.specific_heat=1000.0
rtypes=rtypes+[hotcell]

# define object for atmosphere
top=rocktype('top  ', nad=3, permeability = [perm]*2+[perm],
porosity=poro)
top.conductivity=4 
top.tortuosity=0.0
top.relative_permeability=norp
top.capillarity=nocp
top.specific_heat=1000.0
rtypes=rtypes+[top]

hp=rocktype('hp   ', nad=3, permeability = [perm*10]*2+[perm*10.],
porosity=poro) 
hp.conductivity= 4 #  M/(m K) from Hickey - cotapaxi
hp.tortuosity=0.0
hp.relative_permeability=rp # if single phase this has no effect
hp.capillarity=cp # if single phase this has no effect
hp.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[hp] # add object to list of rocktypes

# define rock types and add cp and rp params
#lp=rocktype('lp   ', nad=3, permeability = [lowk]*2+[lowk],
#porosity=lowporo) 
#lp.conductivity=2 # 3 W/(m K) from Hickey - cotapaxi 
#lp.tortuosity=0.0
#lp.relative_permeability=rp
#lp.capillarity=cp
#lp.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
#rtypes=rtypes+[lp] # add object to list of rocktypes  - not used

main.conductivity=b.conductivity=source.conductivity=top.conductivity=hp.conductivity=conds # reset all conductivities
#%%

if np.size(geo.columnlist) > 1: # can be used to find lateral boundaries in a 2D model - NOTE: WILL NOT WORK FOR 3D
   newlist=np.array([(col,col.centre[0]) for col in geo.columnlist])
   ecol=[newlist[newlist[:,1].argsort()][-1,0]]
   print ecol
else: # if the column list length is only 1 then there can be no lateral boundary.
    ecol=[] # set boundary columns to none

grid=ipt.icegrid(geo,dat,rtypes,ecol,infax=False)#, hpregion=[[0,0,3000],[250,0,6000]])#,heatsource=[[0,0,3000],[1500,0,3050]])
ptg.makeradial(geo,grid,width=width)

## Create TOUGH input file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
#%%
#~ Create sources and sinks GENER block in TOUGH input file.      
#tflux=0.09 # J/m2/s example heat flux


# additional output parameters 
dat.parameter['max_timestep']=3.0e10 # maximum timstep length
dat.parameter['print_interval']=1000 # print (output) frequency to flow.out
dat.parameter['timestep']=[1.0]#[1.0,1000.0] # initial timestep?
dat.output_times['time']=[1.0,3.1558e+08,3.1558e+09,3.1558e+10]#[1.0,1000.0,3.1558e+08,3.1558e+09,3.1558e+10] # predefined output times
dat.output_times['num_times_specified']=len(dat.output_times['time'])
dat.output_times['num_times']=len(dat.output_times['time'])
#dat.parameter['tstop']=1E3*yrsec
dat.output_times['num_times']=20
dat.output_times['time_increment']= 1000*yrsec
#
dat.clear_generators()
ipt.heatgen(mod,geo,dat,grid,heat_flux,function={'type':'log','points':[[1,2.5],[10000,0.24]]})
ptg.gen_constant(mod,geo,grid,dat,constant=1.5e-5,enthalpy=8440.)

geo.write(mod+'/grd.dat')   
# Add data stored in grid object to pyTOUGH dat object - this is what gets turned into the TOUGH input file        
dat.grid=grid
#       
## write vtk of input information - can be view with para view
grid.write_vtk(geo,mod+'/'+mod+'_initial.vtk') 
#   
## write tough2 input file   
dat.write(mod+'/flow2.inp')
        
shutil.copy('dev_files/initial_it2file',mod+'/'+mod)
print time.clock()-t0
