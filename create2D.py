# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 11:30:23 2015
Create Radial model 
@author: glbjch
"""

import pytoughgrav as ptg
import os
from t2grids import *
from t2data import *
import time 
import shutil
from scipy import interpolate

t0=time.clock()
os.chdir("C:\Users\glbjch\Local Documents\Work\Modelling\Gravpaper") # define working directory
mod='20150814_2'
pseudo_topsurf=True
if not os.path.exists(mod):
    os.makedirs(mod)
    
width=10.    
zcells=[2]*50+[10]*20+[2]*55+[10]*4 # 400 m # [10]*34+[2]*80+[10]*19+[2]*30+[10]*5 # 800 m (down to 50 m bsl) # 
surf=ptg.topsurf('dev_files/flatprofile.txt',delim='\t',headerlines=1,width=width)
modelorigin=(0,0)#(586034.886,1852660.465)

## define well locations
#stations=np.genfromtxt('dev_files/Steffie_station_locs.txt', delimiter=',', dtype=None, skiprows=1, usecols=(0,1) ,names='x,y')
#station_dists=[np.sqrt(((modelorigin[0]-x)**2 + (modelorigin[1]-y)**2)) for x,y in zip(stations['x'],stations['y'])]
#station_dists.sort()
#wellx=station_dists
wellx=[50.0,500.0,1100.0,1750,2650,3200]##5.0,50.0,500,1750,2200,2650#######################################################
#wellx=[50.0]################################################################
welly=[width/2.0]*len(wellx) # make y coords same length as x coords
wells=np.hstack((np.transpose([wellx]),np.transpose([welly])))

origin=([0,0,400]) # ([0,0,750])
geo=ptg.geo2D( mod, length=3600, width=width, celldim=10., origin=origin,
              zcells=zcells, surface=surf, wells=wells )


### Create TOUGH input file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
dat=t2data('dev_files/initialflow2.inp') # read from template file 

# define relative permeability and cp paramters to use
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}
norp={'type':5, 'parameters':[]}
cp={'type':11, 'parameters':[0.0,5000.0,-0.001618,0.85,None,None,0.0]}
nocp={'type':1, 'parameters':[0.0,0.0,1.0]}
highk=5.0e-13
lowk=1.0e-16
lowporo=0.1
highporo=0.34

## define rock types - this just generates rock types they are not necessarily assigned to elements - that happens later 
rtypes=[] # creates empty list
dat.incon.clear()
# define object name for rock type. e.g. hp is the python object 'main ' will be the name of the ROCK in TOUGH input file. THE SPACE IN THE NAME IS IMPORTANT - MUST BE 5 char!

# define rock types and add cp and rp params
lp=rocktype('lp   ', nad=3, permeability = [lowk]*2+[lowk],
porosity=lowporo) 
lp.conductivity=1.5 # 3 W/(m K) from Hickey - cotapaxi - not used
lp.tortuosity=0.0
lp.relative_permeability=rp
lp.capillarity=cp
lp.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[lp] # add object to list of rocktypes  - not used

hp=rocktype('hp   ', nad=3, permeability = [highk]*2+[highk],
porosity=highporo) 
hp.conductivity=1.5 # 3 M/(m K) from Hickey - cotapaxi
hp.tortuosity=0.0
hp.relative_permeability=rp # if single phase this has no effect
hp.capillarity=cp # if single phase this has no effect
hp.specific_heat=1000.0 # J/(kg K) from Hickey - cotapaxi
rtypes=rtypes+[hp] # add object to list of rocktypes

# define object for boundary rocktype - not used in current example.....
b=rocktype('nocp ', nad=3, permeability = [highk]*2+[highk],
porosity=highporo)
b.conductivity=1.5 
b.tortuosity=0.0
b.relative_permeability=norp
b.capillarity=nocp
b.specific_heat=1000.0
rtypes=rtypes+[b]

# define object for base rock type - these rocktypes may all be the same but it can be useful to divide them up for defining initial conditions for different regions - there may be other ways to do this.....
#base=rocktype('base ', nad=3, permeability = [highk]*2+[highk],
#porosity=poro)
#base.wet_conductivity=3 
#base.tortuosity=0.0
#base.relative_permeability=norp
#base.capillarity=nocp
#base.specific_heat=1000.0
#rtypes=rtypes+[base]

# define object for top rock type
top=rocktype('atmos', nad=3, density=1.225, permeability = [highk]*2+[highk],
porosity=1.0)
top.conductivity=1.5 
top.tortuosity=0.0
top.relative_permeability=norp
top.capillarity=nocp
top.specific_heat=1000.0
rtypes=rtypes+[top]

lpregion=[[0,0,-50],[1400,0,250]]#[[0,0,500],[1400,0,500]]#
newlist=np.array([(col,col.centre[0]) for col in geo.columnlist])
ecol=[newlist[newlist[:,1].argsort()][-1,0]]
print ecol
#ecol=[geo.columnlist[-1]] # create list of boundary columns (last)

# send to grid2D function to create grid and add rocktype information and define intial conditions 
# also add permeability modifications........
grid = ptg.grid2D(mod,geo,dat,rtypes,ecol,lpregion=lpregion) 

ptg.makeradial(geo,grid,width=width)

dat.parameter['print_block']='ay 40'#############################################################
# add rocktype, element and connection data to dat class
dat.grid=grid


# Define GENER block
#fpms=7.7354e-6 # flux per meter squared
fm=3.64742e-8 # from 1999 - 2013 recharge model within radial model radius. 3.24e-8 # old chapt6 value
fc=2.15803e-6 # from 1999 - 2013 recharge model within radial model radius. -7.199e-7 # old chapt6 value
mingen=2.0e-7 # with fc positive shouldnt be used......

# define constant generation based on elevation relationship.....
#ptg.gen_constant(mod,geo,grid,dat,elev_m=fm,elev_c=fc,mingen=mingen)
if pseudo_topsurf:
    topsurf=np.loadtxt('dev_files/2Dprofile.txt',delimiter='\t',skiprows=1)
    x=topsurf[:,0]
    z=topsurf[:,1]
    s=interpolate.UnivariateSpline(x,z)
    xnew=np.sort([col.centre[0] for col in geo.columnlist])
    znew=s(np.sort(xnew))
    plt.figure()
    plt.plot(x,z,xnew,znew)
    topsurf=np.vstack((xnew,znew)).T
    shutil.copy('dev_files/2Dprofile.txt',mod+'/')

ptg.gen_constant(mod,geo,grid,dat,elev_m=fm,elev_c=fc,mingen=mingen,
                 pseudo_elev=None, pseudo_topsurf=topsurf)
       
## write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk')
#
#geo.write(mod+'/grd.dat')   
## write tough2 input file   
dat.write(mod+'/flow2.inp')
shutil.copy('dev_files/initial_it2file',mod+'/'+mod)
shutil.copy('dev_files/20150527_2_rand.dat',mod+'/')

print time.clock()-t0
   






