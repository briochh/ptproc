# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/briochh/.spyder2/.temp.py
"""

from t2data import *
from t2grids import *
from t2listing import *
import os
import pytoughgrav as ptg
import matplotlib.pyplot as plt
import matplotlib.colors as mcols

plt.close('all')
os.chdir('/Users/briochh/Documents/Workhere/testing')
mod='topoRAD_4_0__0_24surfperm_1E-14_FLATTOP' # define model name
os.chdir(mod)
if not os.path.exists('results'): 
    os.makedirs('results')   
    
dat=t2data('flow2.inp')
grid=dat.grid
geo=mulgrid('grd.dat')
ptg.makeradial(geo,None,1.0)
#results=t2listing('flow2.out')   # comment out if results has already been loaded........
os.chdir('results')
if not os.path.exists('surfaceflow_test'):
    os.makedirs('surfaceflow_test')

# find atmosphere blocks
grid2=t2grid().fromgeo(geo) # grid saved in flow2.inp does not contain atmos information required.
atmos=grid2.atmosphere_blocks
# find atmosphere connections
atmosconn=[cell.connection_name for cell in atmos] # returns a list of sets - !!!! sets are not ordered so cannot be indexed !!! 

flowH={} # if was we to record flow by connection in a dictionary 
flowF={}
totqh=glac_toqh=totqf=glac_toqf=np.zeros(results.num_times) # arrays for adding total flows at each timstep 

X=np.array([]).reshape(0,1) # array for storing the x dimension of each connection 
qhts=[] # empty list for storing all flows
qfts=[]
for atmosset in atmosconn: # loop over each set of atmos connections
    for conn in atmosset: # loop over each connection in set
        cell1,cell2=conn
        [(tq,qh),(tqf,qf)] = results.history([('c', conn, 'FLOH'),('c',conn,'FLOF')]) # here just pull out flow of HEAT (could do more)
        # flow from cell1 to cell2 is negative we want to make it uniform so that flow into the atmosphere is always negative
        if grid2.block[cell2] in atmos:
            qh=qh/grid2.connection[conn].area
            qf=qf/grid2.connection[conn].area
            cen=grid2.block[cell1].centre # 
        elif grid2.block[cell1] in atmos:
            qh=-qh/grid2.connection[conn].area # if atmos cell is cell1 the sign of the flow needs to be switched
            qf=-qf/grid2.connection[conn].area            
            cen=grid2.block[cell2].centre
        flowH[conn]=cen[0],qh # dictionary of flow timseries for each connection
        flowF[conn]=cen[0],qf        
        # compile all the flow time series (for each connection)        
        if qhts == []:
            qhts=[qh] # the first time through need to initialise our list
        else:
            qhts=np.concatenate((qhts,[qh]),axis=0) # add the current flow time series to our compilations
        if qfts == []:
            qfts=[qh] # the first time through need to initialise our list
        else:
            qfts=np.concatenate((qfts,[qf]),axis=0) # add the current flow time series to our compilations
        
        X=np.append(X,cen[0]) # also create array of connection distannce
        totqh=np.add(totqh,qh) # timeseries of total heat flux into atmos 
        totqf=np.add(totqf,qf) # timeseries of total fluid flow into atmos 
        # can look at total heat flux that is within a certain X limit (e.g. within the range of a glacier)        
        if cen[0] <= 2000.0:
            glac_toqf=np.add(glac_toqf,qf) 
            glac_toqf=np.add(glac_toqf,qf) 
        #atxs[conn]=cen[0]

pcmap = 'seismic'
minimum_log_level = 1e-9
maximum_scale_level = 4e-5
# Create a 'logarithmic' data normalization.
pnorm = mcols.SymLogNorm(linthresh=minimum_log_level,
                                 linscale=0,
                                 vmin=-maximum_scale_level,
                                 vmax=maximum_scale_level)
# a quick plot of flows into atmosphere at X and time.
plt.figure()
plt.pcolormesh(X,tq,qhts.T, vmin=-0.5, vmax=0.5)
plt.colorbar().set_label(r'Flow of heat into atmosphere 'r'(W/m$^{2}$)')
plt.xlim((X.min(),X.max()))
plt.ylim((1e10,tq.max()))
plt.show()        

plt.figure()
plt.pcolormesh(X,tqf,qfts.T, cmap=pcmap, norm=pnorm)
plt.colorbar(format='%.0e').set_label(r'Flow of fluid into atmosphere 'r'(kg/s/m$^{2}$)')
plt.xlim((X.min(),X.max()))
plt.ylim((1e10,tqf.max()))
plt.show()   

plt.figure()
plt.pcolormesh(X,tqf,qfts.T, vmin=-4e-5, vmax=4e-5)
plt.colorbar(format='%.0e').set_label(r'Flow of fluid into atmosphere 'r'(kg/s/m$^{2}$)')
plt.xlim((X.min(),X.max()))
plt.ylim((1e10,tqf.max()))
plt.show()      

plt.figure()
plt.plot(X,qfts[:,-1])
plt.yscale('symlog', linthreshy=minimum_log_level)

#results.write_vtk(geo,'output.vtk',grid=grid,lows=True)
