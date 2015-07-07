





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
import scipy as spy

def icegrid(geo,dat,rocks,boundcol,lpregion=None,hpregion=None,heatsource=None,satelev=0.0,atmosP=1.013e5,pmx_lamda=0.004, glacier_limit=2500., infax=False):
    """
    Method for defining ice grid. Varies slightly from ptg method.
    """
    grid = t2grid().fromgeo(geo)
    atmos=grid.atmosphere_blocks   
    # loop over list of rocktypes and add to pytough gird object
    for rock in rocks:
        grid.add_rocktype(rock) 
    # loop over every element block in the grid
    for blk in grid.blocklist[0:]:
        blk.hotcell=False
        lay=geo.layer[geo.layer_name(str(blk))] # layer containing current block
        col=geo.column[geo.column_name(str(blk))] # column containing current block
        tlay=geo.column_surface_layer(col)    
        hmax=geo.block_surface(tlay,col)
        if blk in atmos:
            rocktype='top  ' # assign rocktype "top  "
            initP=atmosP*spy.power(1.-(hmax*2.25577e-5),5.25588)  # initial presure condition - MAY NOT BE APPROPRIATE - WHAT IS THE PRESSURE UNDER THICK GLACIER AT 5000 m amsl??         
            Tmin=2.      
            if blk.centre[0] <= glacier_limit: 
                initT=Tmin # initial temperature - TOUGH2 doesn't seem to like < 1.0 C
            else:
                initT = 25.8 - (hmax*(5.4/1000.)) # 15.+((2000.-blk.centre[2])*(5.4/1000.0))
                if initT <= Tmin: initT=Tmin
            if (hpregion is not None and 'hp   ' in grid.rocktype.keys()):
                for hpr in hpregion.values():
                    if (blk.centre[2] > hpr[0][2] and 
                    blk.centre[2] <= hpr[1][2] and 
                    blk.centre[0] > hpr[0][0] and 
                    blk.centre[0] <= hpr[1][0]): #if in hp region
                        rocktype='hp   ' # this allows a different pmx for atmos above highperm
            initSG=0.999 # initial gas saturation
            infvol=False # already given 1e50 volume
            pmx=grid.rocktype[rocktype].permeability[0]
            rocktype='top  ' # resets to rocktype "top  "
        else:
            rocktype = 'main '
            initP=atmosP*spy.power(1.-(hmax*2.25577e-5),5.25588)+(765.*9.81*abs(hmax-blk.centre[2]))
            initSG=0.0
            initT=Tmin+((np.abs(hmax-blk.centre[2])/100.0)*3.0)
            infvol=False
            if lay==geo.layerlist[-1]:
                rocktype='sourc'
            if (lpregion is not None and 'lp   ' in grid.rocktype.keys() and
                blk.centre[2] > lpregion[0][2] and 
                blk.centre[2] <= lpregion[1][2] and 
                blk.centre[0] > lpregion[0][0] and 
                blk.centre[0] <= lpregion[1][0]): # if in lp region
                rocktype='lp   '     
            if (hpregion is not None and 'hp   ' in grid.rocktype.keys()):
                for hpr in hpregion.values():
                    if (blk.centre[2] > hpr[0][2] and 
                    blk.centre[2] <= hpr[1][2] and 
                    blk.centre[0] > hpr[0][0] and 
                    blk.centre[0] <= hpr[1][0]): #if in hp region
                        rocktype='hp   '
            if (heatsource is not None and 
                blk.centre[2] > heatsource[0][2] and 
                blk.centre[2] <= heatsource[1][2] and 
                blk.centre[0] > heatsource[0][0] and 
                blk.centre[0] <= heatsource[1][0]): # if in heatsource region
                rocktype='hotcl'
                initT=350
                infvol=True
                blk.hotcell=True
            if infax is True and col is geo.columnlist[5] and lay == tlay:
                print "inf vol top axis cell " + blk.name
                infvol=True
                initSG=10.9999
                rocktype='bound'
            if col in ecol:
                rocktype='bound'
                infvol=True
                initP=atmosP*spy.power(1.-(hmax*2.25577e-5),5.25588)+(767.*9.81*abs(satelev-blk.centre[2]))
                if blk.centre[2]>satelev:
                    initSG=0.999
                    initP=atmosP*spy.power(1.-(hmax*2.25577e-5),5.25588)
            pmx=ipt.pmxcalc(blk,grid,hmax,rocktype,0.004,800.)      
        ptg.rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx,infvol=infvol)
    return grid

#%% Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0=time.clock()
os.chdir("C:\Users\glbjch\Local Documents\Work\Modelling\Cotapaxi") # define working directory
mod='Boiltest_20150706_3'
print mod
if not os.path.exists(mod):
    os.makedirs(mod)
#%%
yrsec=3600*365.25*24
origin=[0,0,5500] # position of first cell in space
width=1.0
zcells=[10]*30+[50]*4+[100]*10+[50]*6+[25]*6+[10]*5    #+[100]*6+[50]*10+[10]*20    
dy=1 # size of cell in y direction
xcells=[5]*12+[10]*54  #+[100]*6+[50]*10+[10]*20

#surf=ptg.topsurf('dev_files/Topography_crater_f.txt',delim='\t',headerlines=1,width=width)

geo=ipt.icegeo( mod, width=width, celldim=10., origin=origin,
              zcells=zcells, xcells=xcells, atmos_type=1, min_thick=2.0)

#%%
dat=t2data('dev_files/initialflow2.inp') # read from template file 
dat.parameter['print_block']=' w 46' # define element to print in output - useful for loggin progress of TOUGH sim 
dat.multi['num_equations']=3 # 3 defines non isothermal simulation

perm=5.0e-13 # define permeability
poro=0.1  # define porosity
rp={'type':1, 'parameters':[0.3,0.1,0.9,0.7]}#{'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]} # relative permeability functions and parameters - if single phase not necessary  
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
source.relative_permeability=rp
source.capillarity=cp
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

hp=rocktype('hp   ', nad=3, permeability = [perm]*2+[perm*10.],
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

grid=icegrid(geo,dat,rtypes,ecol,satelev=5300.)#, hpregion={'hpr1':[[0,0,3000],[50,0,6000]],'hpr2':[[200,0,3000],[250,0,6000]]})#,heatsource=[[0,0,3000],[1500,0,3050]])
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
ipt.heatgen(mod,geo,dat,grid,heat_flux,function={'type':'log','points':[[1.,5.],[10000.,0.24]]})
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
