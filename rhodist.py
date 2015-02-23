# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\glbjch\.spyder2\.temp.py
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os
import time
import scipy.constants
import scipy.io as io

plt.close('all')

save='no'
savevtk='no'
t0=time.clock()
mod='20140827_2_py_it'


os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/'+mod)

# read input grid and data
geo=mulgrid('2dgrd.dat')
dat=t2data('flow2.inp')
grid=dat.grid # define input grid

# define well locations
#wellx=[50,500,1750,2200,2650]
wellx=[0.5]
welly=[0.5]*len(wellx)
wells=hstack((np.transpose([wellx]),np.transpose([welly])))

t1=time.clock()
t=t1-t0
print 'time2setup=',t

# read output file
t0=time.clock()
#results=t2listing('flow2.out') # read output file
t1=time.clock()
t=t1-t0
print 'time2read .out=',t

t0=time.clock()
density=997.0479
cen=np.array([blk.centre for blk in grid.blocklist[1:]])
dz=np.array([geo.layer[geo.layer_name(blk.name)].thickness for blk in grid.blocklist[1:]])
phi=np.array([blk.rocktype.porosity for blk in grid.blocklist[1:]])
vol=np.array([blk.volume for blk in grid.blocklist[1:]])
for i in range(vol.size):
   if vol[i]>1e20:
       vol[i]=vol[i-1]
   
cenvoldz={'xyz':cen,'volume':vol,'thickness':dz}
t1=time.clock()
t=t1-t0
print 'time to produce cell geom info',t


t0=time.clock()
results.first()       
data={}
count=0
for i in results.times:
   print('timestep %d out of %d' % (count+1,results.times.size))
#   t=results.time 
   sat=results.element['SL'][1:]
   drho=phi*sat*density
   #   times=np.empty(drho.shape)
   #   times.fill(time)
   data['t'+str(count)]=drho
#   ts['t'+str(count)]=t
   results.next()
   count=count+1
   
t1=time.clock()
t=t1-t0
print 'time to calculate density changes for every ts=',t
    
t0=time.clock()
io.savemat('rho',data)
io.savemat('times',dict(times=results.times))
io.savemat('cenvoldz',cenvoldz)
t1=time.clock()
t=t1-t0
print 'time2save .mat files',t
#np.savetxt('CenRhoVolZ.dat',np.concatenate((cen,np.reshape(drho,[drho.shape[0],1]),np.reshape(vol,[vol.shape[0],1]),np.reshape(dz,[dz.shape[0],1])),1))


#results.times.shape
#if not os.path.exists('results'):
#    os.makedirs('results')
#os.chdir('results')
#

#yrsec=3600*24*365.35
#
#wellno=0
#for well in wells:
#    wellno=wellno+1
#    print well
#    t0=time.clock()
#    col=geo.column_containing_point(well)
#    wellblk=[]
#    zlist=[]
#    for lay in geo.layerlist[2:]:
#        if geo.block_name(lay.name,col.name) in geo.block_name_list:
##            print lay.name
#            blk=grid.block[geo.block_name(lay.name,col.name)]
##            #blk.sl=results.history(('e',blk.name,'SL'))
##            wellsl=np.concatenate((wellsl,[results.history(('e',blk.name,'SL'))
##                                  [1]]),axis=0)
##            
#            wellblk.append(blk)
#    t1=time.clock()
#    t=t1-t0
#    print 'time4wellblk',t 
#    
#    t0=time.clock()
#    selsl=[]
#    outsl=[]
#    
#    outd={}
#    for blk in wellblk:
##        print blk
#        selsl=selsl+[('e',blk.name,'SL')]
#        zlist=np.concatenate((zlist,[blk.centre[2]]))
#    outsl=results.history(selsl) # saturation of all cells in column over time
#    times=outsl[0][0]
#    wellsl=np.array([]).reshape(0,len(times))
#    wellro=np.array([]).reshape(0,len(times))
#    well_water_mass=np.zeros(len(times))
#    wellvol=0
#    
#    for blk,outsl1 in zip(wellblk,outsl):
#        outd[blk.name]=outsl1[1]
#        wellsl=np.concatenate((wellsl,[outsl1[1]]))
#        wellro=np.concatenate((wellro,[outsl1[1]*density*blk.rocktype.porosity]))
#        if blk is not wellblk[-1]:
#           wellvol=wellvol+blk.volume
#           well_water_mass=well_water_mass+(outsl1[1]*blk.rocktype.porosity*blk.volume*density)
#        #well_water_mass=well_water_mass+(outsl1[1]*blk.rocktype.porosity*blk.volume*density)
#        
#    microgal=well_water_mass*2*pi*scipy.constants.G/col.area*10**8   
#    microgal=microgal-microgal[0]
#    t1=time.clock()
#    t=t1-t0
#    print 'time4well',t
#     
#    t0=time.clock()
#    im1=plt.figure()
#    plt.plot(times/3.15576e7,microgal)
#    plt.ylabel(r'$\Delta g$ (microgal)')
#    plt.xlabel('Time (years)')
#    plt.axis([0.0, times.max()/yrsec,None,None])
#    
#    T,Z=np.meshgrid(times/3.15576e7,zlist)    
# #   plt.figure()
#  #  im=plt.imshow(wellsl, interpolation='bilinear',origin='upper',
#  #                vmin=wellsl.min(),vmax=wellsl.max(), cmap=cm.jet_r,
#   #               extent=(0.0,6.2221e+09,-250,750),aspect='auto')
#                  
#    im2=plt.figure()
#    profplt=plt.pcolormesh(T,Z,wellsl,cmap=cm.jet_r,vmin=0.0,vmax=1.0,shading='flat')
#    plt.axis([T.min(), T.max(), Z.min(), Z.max()])
#    plt.ylabel('Z (m)')
#    plt.xlabel('Time (years)')    
#    cbar=plt.colorbar(profplt,orientation='vertical')
#    cbar.set_label('Saturation')   
#    t1=time.clock()
#    t=t1-t0
#    print 'time2plotwell',t 
#
#
#    if save is 'yes':  
#       t0=time.clock()      
#       zt_density_matrix=np.concatenate((
#       [np.concatenate((np.array([0]),times))],
#        np.concatenate((zlist.reshape(len(zlist),1),wellro),
#                       axis=1)),
#                       axis=0)
#       f = open('resultxt.txt','w')
#       f.write('Model = '+mod+'\n'
#               'Mass max (kg) =' +str(well_water_mass.max())+'\n'
#               'Mass min (kg) =' +str(well_water_mass.min())+'\n'
#               'Max amplidute (mass)='+str(well_water_mass.max()-well_water_mass.min())+'\n'
#               'grav max (microgal) =' +str(microgal.max())+'\n'
#               'grav min (microgal) =' +str(microgal.min())+'\n'
#               'Max amplidute (grav)='+str(microgal.max()-microgal.min())+'\n')
#       f.close()
#       savetxt('ztro.dat',zt_density_matrix)                
#       savetxt('waterweight'+str(wellno)+'.dat',zip(times,well_water_mass))
#       savetxt('microgal'+str(wellno)+'.dat',zip(times,microgal))
#       im2.savefig('sl_t_profile'+str(wellno)+'.png',dpi=300)
#       im2.savefig('sl_t_profile'+str(wellno)+'.eps')
#       im1.savefig('microgal'+str(wellno)+'.pdf')
#
#    #savefig('sl_t_profile.pdf')
#       t1=time.clock()
#       t=t1-t0
#       print 'time2saveplot',t
#       
#       
##for i in range(0,len(times)):
# #   t=times[i]
#  #  profilez=[]
#   # profilesl=[]
#    #slgrd=np.concatenate((blk.sl[1]) for blk in wellblk[1:])
#     #  profilez.append(blk.centre[2])
#      # profilesl.append(blk.sl[1][i])
#       #plot(profilesl,profilez)
#       
#      
##thingys=[wellblk[1].sl[1]]
##for thing in wellblk[2:]:
# # thingys=np.concatenate((thingys,[thing.sl[1]]))      
#        
#if savevtk is 'yes':
#   t0=time.clock()
#   results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
#   t1=time.clock()
#   t=t1-t0
#   print 'time2writevtks',t
