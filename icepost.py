# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 17:48:26 2015

@author: glbjch
"""

from t2data import *
from t2grids import *
from t2listing import *
import os
import pytoughgrav as ptg
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Cotapaxi')
mod='Cota20150604_1' # define model name
os.chdir(mod)
if not os.path.exists('results'): 
    os.makedirs('results')   
    
dat=t2data('flow2.inp')
grid=dat.grid
geo=mulgrid('grd.dat')
ptg.makeradial(geo,None,1.0)
results=t2listing('flow2.out')   # comment out if results has already been loaded........
#os.chdir('/Users/molly/Documents/Project/topoRAD/Results_infvol_1000')
#if not os.path.exists('surfaceflow_test'):
#    os.makedirs('surfaceflow_test')

# find atmosphere blocks
grid2=t2grid().fromgeo(geo) # grid saved in flow2.inp does not contain atmos information required.
atmos=grid2.atmosphere_blocks
# find atmosphere connections
atmosconn=[cell.connection_name for cell in atmos] # returns a list of sets - !!!! sets are not ordered so cannot be indexed !!! 

flowH={} # if was we to record flow by connection in a dictionary 
totq=glac_toq=np.zeros(results.num_times) # arrays for adding total flows at each timstep 

X=np.array([]).reshape(0,1) # array for storing the x dimension of each connection 
Area=np.array([]).reshape(0,1)
qts=[] # empty list for storing all flows
for atmosset in atmosconn: # loop over each set of atmos connections
    for conn in atmosset: # loop over each connection in set
        cell1,cell2=conn
        tq,q = results.history([('c', conn, 'FLOH')]) # here just pull out flow of HEAT (could do more)
        # flow from cell1 to cell2 is negative we want to make it uniform so that flow into the atmosphere is always negative
        if grid2.block[cell2] in atmos:
            q=q/(grid2.connection[conn].area)
            cen=grid2.block[cell1].centre # 
        elif grid2.block[cell1] in atmos:
            q=-q/(grid2.connection[conn].area) # if atmos cell is cell1 the sign of the flow needs to be switched
            cen=grid2.block[cell2].centre
        flowH[conn]=cen[0],q # dictionary of flow timseries for each connection
        # compile all the flow time series (for each connection)        
        if qts == []:
            qts=[q] # the first time through need to initialise our list
        else:
            qts=np.concatenate((qts,[q]),axis=0) # add the current flow time series to our compilations
        X=np.append(X,cen[0]) # also create array of connection distannce
        Area= np.append(Area,grid2.connection[conn].area)
        totq=np.add(totq,q) # timeseries of total heat flux into atmos 
        # can look at total heat flux that is within a certain X limit (e.g. within the range of a glacier)        
        if cen[0] <= 2500.0:
            glac_toq=np.add(glac_toq,q) 
        #atxs[conn]=cen[0]










## a quick plot of flows into atmosphere at X and time.
#plt.figure()
#plt.pcolormesh(X,tq/3600/24/365.25,qts.T, rasterized=True)
#plt.colorbar().set_label(r'Flow of heat into the surface 'r'(W/m$^{2}$)')
##plt.xlim((0,2500)
##plt.ylim((0 , (tq/3600/24/365.25).max()))
##plt.title('Flow of heat into the glacier following perturbation')
##plt.xlabel('
##plt.ylabel('
#
#plt.savefig("/Users/molly/Documents/Project/Restart_images/Post_heating/" + mod + "Heat flow.pdf", dpi=500)
#plt.close()    
#
#
#
##results.write_vtk(geo,'output.vtk',grid=grid,lows=True)
meltratematrix= (qts.T/-3.35E5)
#plt.figure()
#plt.pcolormesh(X,tq/3600/24/365.25,meltratematrix, rasterized=True)
#plt.colorbar().set_label(r'Melt rate 'r'(Kg/y/m$^{2}$)')
##plt.xlim((0,2500))
##plt.ylim((0 , (tq/3600/24/365.25).max()))
##plt.title('Glacial melt rate following perturbation')
##plt.xlabel('
##plt.ylabel('
#
#plt.savefig("/Users/molly/Documents/Project/Restart_images/Post_heating/" + mod + "Melt rate.pdf", dpi=500)
#plt.close()


#############

i=0
meltrate=np.zeros(len(tq))
for t in tq:
    for x,A,r in zip(X,Area,meltratematrix[i]):
        if r> 0 and x < 2500:
            meltrate[i]= meltrate[i] + (r*A)
    i=i+1

meltrate_mpyr= (meltrate*(3600*24*365.25)/1000)/(np.pi*(2500^2))

#np.savetxt("melt_time.txt", np.vstack(([tq], [meltrate],[meltrate_mpyr])).T)
   
#plt.figure()
#plt.plot(tq/3600/24/365.25,meltrate_mpyr)
#plt.xlim(0, 30000)
#plt.xlabel('Years after end of thermal perturbation')
#plt.ylabel('Average melt rate (m/yr) at glacial base')
#plt.title('Average basal meltrate following thermal perturbation \n of 350 degrees for 1ky at 2.5km depth')
#plt.savefig("/Users/molly/Documents/Project/Restart_images/Meltrate_following_perturbation*", dpi=500)
#plt.close() 


#==============================================================================
# Change in meltrate over time at various X
#==============================================================================

meltratechange=np.array([tstep-meltratematrix[0] for tstep in meltratematrix])
plt.figure()
plt.pcolormesh(X,np.log10(tq),qts.T, rasterized=True, cmap='rainbow')
#plt.xlim((0,2500))
#plt.ylim((0 , 30000))
plt.xlabel('Distance from centre of volcano (m)')
plt.ylabel('Years after end of thermal perturbation')
cbar=plt.colorbar()
cbar.set_label(r'Increase in melt rate 'r'(Kg/y/m$^{2}$)')
plt.title('Change in meltrate following thermal perturbation \n with distance from crater centre')
#plt.savefig("/Users/molly/Documents/Project/Restart_images/Meltrate_at_X_with_time*", dpi=500)
#plt.close()