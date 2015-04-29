# -*- coding: utf-8 -*-
"""Reading iTOUGH results"""
"""
Created on Thu May 22 09:00:37 2014

@author: glbjch
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.tri as tri
import os
import time
import scipy.constants
import scipy.integrate as integrate
import scipy.interpolate

plt.close('all')

read='yes'#################################################################
save='yes'##################################################################
savevtk='yes'###############################################################
t0=time.clock()
mod='20140905_py_it'


os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/'+mod)

## read input grid and data
if read is 'yes':
   geo=mulgrid('2dgrd.dat')
   dat=t2data('flow2.inp')
   grid=dat.grid # define input grid

## define well locations
wellx=[50.0,500,1750,2200,2650]##5.0,50.0,500,1750,2200,2650#######################################################
#wellx=[0.5]################################################################
welly=[0.5]*len(wellx) # make y coords same length as x coords
wells=hstack((np.transpose([wellx]),np.transpose([welly])))
  
t1=time.clock()
t=t1-t0
print 'time2setup=',t

## read output file
if read is 'yes':
   t0=time.clock()
   results=t2listing('flow2.out') # read output file
   t1=time.clock()
   t=t1-t0
   print 'time2read .out=',t

## create directory for results
if not os.path.exists('results'):
    os.makedirs('results')
os.chdir('results')

density=997.0479 # density of water
yrsec=3600*24*365.35 # 1 year in seconds


## create arrays conatining information about the geometry of each block. 
## this is quicker than acessing the listings class every time
dz=dx=[]
cen=np.array([]).reshape(0,3)
for blk in grid.blocklist[1:]:
    dz=np.concatenate((dz,[geo.layer[geo.layer_name(blk.name)].thickness])) # array of thickness of each block
    dx=np.concatenate((dx,[geo.column[geo.column_name(blk.name)].side_lengths[1]])) # array of x direction cell widths
    cen=np.concatenate((cen,[blk.centre]))

xzarea=dz*dx # array of cell x by z area
xs=np.array([c[0] for c in cen])
zs=np.array([c[2] for c in cen])
#X,Z=np.meshgri
## numerical solution of ring integral from 0 to 2pi slightly quicker
## set up arrays of theta increments
ntheta=1000.0 # number of increments
dtheta=2.0*pi/ntheta # length of element of integral
thetas=np.arange(0.0,2.0*pi,dtheta) #array of integral element lengths


## loop over each gravity survey point given in wellx
wellno=0
for well in wells:
    twell=time.clock()
    
    wellno=wellno+1
    print well
    wellblk=[] # list for block names in well column
    zlist=[] # array of cell depths in column
    
    
    ## generate list of cells beneath survey point (within well)
    ## usedfor bouguer slab
    t0=time.clock()
    col=geo.column_containing_point(well) # column that contains the survey point
    for lay in geo.layerlist[1:]:
        if geo.block_name(lay.name,col.name) in geo.block_name_list:
#            print lay.name
            blk=grid.block[geo.block_name(lay.name,col.name)]
#            #blk.sl=results.history(('e',blk.name,'SL'))
#            wellsl=np.concatenate((wellsl,[results.history(('e',blk.name,'SL'))
#                                  [1]]),axis=0)
#            
            wellblk.append(blk) # list of cells in well column (used in bouguer slab calculation)
            #outd[blk.name]=[] # set up dictionary to hold well saturation
            zlist=np.concatenate((zlist,[blk.centre[2]]))
    
    t1=time.clock()
    t=t1-t0
    print 'time4wellblk',t
    
    zp=col.surface # surface eleaviton  - denotes z elevation for gravity measurement
    results.first() # find first timestep for listings file
    gravt=[] # set up gravity over time list
    well_water_mass=[] # set up array for mass of water in current survey column (well), over time
    wellsl=np.array([]).reshape(len(wellblk),0) # array for saturation in well over time
    wellro=np.array([]).reshape(len(wellblk),0) # array for water density in well over time
    
    
    ## compute elemental contribution for each element in domain at current gravity survey location
    i=0
    comp=[]
    for blk in grid.blocklist[1:]:
        alpha=blk.centre[0] # radial distance to element
        blkz=blk.centre[2] # elevation of element
        blkxzarea=xzarea[i] # area of element
        # dtheta ingetral leght method:
        # distance to every element from survey point
        #      s=((alpha**2)-(2*alpha*well[0]*cos(thetas))+(well[0]**2)+((zp-blkz)**2))**(1./2.) 
        #       val=(sum(dtheta/s**3)) # fractional contribution of current ring
        # python integrate method: (slower)?
        I= lambda theta: 1/(((alpha**2)-(2*alpha*well[0]*np.cos(theta))+(well[0]**2)+((zp-blkz)**2))**(1.5)) 
        val,err= integrate.quad(I, 0.0, 2*np.pi)
        # multiply by radius of ring, vertical distance from survey point to ring, xzarea 
        comp.append((zp-blkz)*alpha*val*blkxzarea)
        i+=1
        
    
    # loop over results table untill desired time (in years)
    count=0
    while results.time/yrsec <= 110:
        t0=time.clock()
        print('timestep %d out of %d' % (count+1,results.times.size))
        print results.time/yrsec
        sat=results.element['SL'][1:] # pull out saturation index [0] is the atmosphere cell so no use to us  
        dg=[] # array for collecting together elemental contributions of integral method
        twellsl=[] # array for collecting saturation in well for current timestep
        twellrho=[] # array for collecting saturation density in well for current timestep
        twell_water_mass=0 # starting point for water mass in well at current timestep
        # loop over every block for current timestep
        i=0
        for blk in grid.blocklist[1:]: # dont bother with atmos cell
            blkrho=sat[i]*blk.rocktype.porosity*density # saturation density in element
            dg.append(comp[i]*blkrho) # gravity conribution of element 
            i+=1
            ## For bouguer slab.....
            #calculated saturation, saturation density and water mass in column at current timestep
            if blk in wellblk:
                #outd[blk.name].append(sat[i])
                twellsl.append([sat[i]])
                twellrho.append([blkrho])
                if blk is not wellblk[-1]:
                    #wellvol=wellvol+blk.volume
                    twell_water_mass=twell_water_mass+(blkrho*blk.volume)
#                wellsl=np.concatenate((wellsl,[outsl1[1]]))
#                wellro=np.concatenate((wellro,[outsl1[1]*density*blk.rocktype.porosity]))
            i+=1 # inrement to next element
        #compile array of saturation, density and water mass in current column for all timesteps
        wellsl=np.concatenate((wellsl,twellsl),axis=1)
        wellro=np.concatenate((wellro,twellrho),axis=1)
        well_water_mass.append(twell_water_mass)
        
        # Total axissummetric intgral gravity due to water mass in model for current time step
        grav=6.67e-3*sum(dg)
        gravt.append(grav) # compile gravity due to water at current survey point for each timestep 
        print grav
        results.next() # move to next timestep 
        count+=1 # incrment timestep counter
        dt=time.clock()-t0
        print('time for timestep = %e s' % dt)

    # Bouguer slab gravity approximations
    well_water_mass=np.array(well_water_mass)
    microgal=well_water_mass*2*pi*scipy.constants.G*(10**8)/col.area   
    microgal=microgal-microgal[0] # gravity difference    
    
    t1=time.clock()
    t=t1-twell
    print 'time4well_calculations',t    
    
    ## Plot results for current well
    t0=time.clock()
    times=results.times[0:count]/yrsec # convert times calculted to yrs 
    
    # test plot of contibutions
    im=plt.figure(figsize=[8,3.6])
    plt.scatter(xs, zs, c=np.array(dg/xzarea)*6.67e-3, edgecolor='none', marker='s')
    plt.colorbar()
    plt.xlim((xs.min(),xs.max()))
    plt.ylim((zs.min(),zs.max()))
    plt.show()
    
    # integral gravity time series
    intgravplt=figure()
    plt.plot(results.times[0:count]/yrsec,gravt-gravt[0])
    plt.ylabel(r'$\Delta g$ (microgal)')
    plt.xlabel('Time (years)')
    plt.axis([0.0, times.max(),None,None])
    
    # Bouguer slab estimation fomr column saturation
    im1=plt.figure()
    plt.plot(times,microgal)
    plt.ylabel(r'$\Delta g$ (microgal)')
    plt.xlabel('Time (years)')
    plt.axis([0.0, times.max(),None,None])
    
    # column satuation over time
    T,Z=np.meshgrid(times,zlist)
    im2=plt.figure()
    profplt=plt.pcolormesh(T,Z,wellsl,cmap=cm.jet_r,vmin=0.0,vmax=1.0,shading='flat')
    plt.axis([T.min(), T.max(), 0, Z.max()])
    plt.ylabel('Z (m)')
    plt.xlabel('Time (years)')    
    cbar=plt.colorbar(profplt,orientation='vertical')
    cbar.set_label('Saturation')   
    t1=time.clock()
    t=t1-t0
    print 'time2plotwell',t
    
    if save is 'yes':  
       t0=time.clock()      
       zt_density_matrix=np.concatenate((
       [np.concatenate((np.array([0]),times))],
        np.concatenate((zlist.reshape(len(zlist),1),wellro),
                       axis=1)),
                       axis=0)
       f = open('resultxt_'+str(wellno)+'test.txt','w')
       f.write('Model = '+mod+'\n'
               'Mass in col max (kg) =' +str(well_water_mass.max())+'\n'
               'Mass in col min (kg) =' +str(well_water_mass.min())+'\n'
               'Max col amplidute (mass)='+str(well_water_mass.max()-well_water_mass.min())+'\n'
               'grav (Boug) max (microgal) =' +str(microgal.max())+'\n'
               'grav (Boug) min (microgal) =' +str(microgal.min())+'\n'
               'Max (Boug) amplidute (grav)='+str(microgal.max()-microgal.min())+'\n'
               'grav (int_axsym) max (microgal)='+str((gravt-gravt[0]).max())+'\n'
               'grav (int_axsym) min (microgal)='+str((gravt-gravt[0]).min())+'\n'
               'Max (int_axsym) amplitude (microgal)='+str((gravt-gravt[0]).max()-(gravt-gravt[0]).min())+'\n')
       f.close()
       savetxt('ztro.dat',zt_density_matrix)                
       savetxt('waterweight'+str(wellno)+'.dat',zip(times,well_water_mass))
       savetxt('microgal'+str(wellno)+'.dat',zip(times,microgal))
       savetxt('axsym_int_microgal'+str(wellno)+'.dat',zip(times,(gravt-gravt[0])))
       im.savefig('elemental_cont'+str(wellno)+'.pdf')       
       im2.savefig('sl_t_profile'+str(wellno)+'test.png',dpi=300)
       #im2.savefig('sl_t_profile'+str(wellno)+'.eps')
       im1.savefig('microgal'+str(wellno)+'.pdf')
       intgravplt.savefig('axsym_int_grav'+str(wellno)+'.pdf')
#
     #savefig('sl_t_profile.pdf')
       t1=time.clock()
       t=t1-t0
       print 'time2saveplot',t
    
    t1=time.clock()
    t=t1-twell
    print 'total time for well',t
    
#       plt.imshow(np.reshape(dg,(geo.get_num_layers()-2,geo.get_num_columns())))
##plt.gca().invert_yaxis()
#       plt.colorbar()
        
 
    
#    t0=time.clock()
#    selsl=[]
#    outsl=[]
#    
#    
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
#    
#    
#     
#    t0=time.clock()
#    im1=plt.figure()
#    plt.plot(times/yrsec,microgal)
#    plt.ylabel(r'$\Delta g$ (microgal)')
#    plt.xlabel('Time (years)')
#    plt.axis([0.0, times.max()/yrsec,None,None])
#    
#    T,Z=np.meshgrid(times/yrsec,zlist)    
# #   plt.figure()
#  #  im=plt.imshow(wellsl, interpolation='bilinear',origin='upper',
#  #                vmin=wellsl.min(),vmax=wellsl.max(), cmap=cm.jet_r,
#   #               extent=(0.0,6.2221e+09,-250,750),aspect='auto')
#                  
#    im2=plt.figure()
#    profplt=plt.pcolormesh(T,Z,wellsl,cmap=cm.jet_r,vmin=0.0,vmax=1.0,shading='flat')
#    plt.axis([T.min(), T.max(), 0, Z.max()])
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
#       #im2.savefig('sl_t_profile'+str(wellno)+'.eps')
#       im1.savefig('microgal'+str(wellno)+'.pdf')
#
#    #savefig('sl_t_profile.pdf')
#       t1=time.clock()
#       t=t1-t0
#       print 'time2saveplot',t
       
       
#for i in range(0,len(times)):
 #   t=times[i]
  #  profilez=[]
   # profilesl=[]
    #slgrd=np.concatenate((blk.sl[1]) for blk in wellblk[1:])
     #  profilez.append(blk.centre[2])
      # profilesl.append(blk.sl[1][i])
       #plot(profilesl,profilez)
       
      
#thingys=[wellblk[1].sl[1]]
#for thing in wellblk[2:]:
 # thingys=np.concatenate((thingys,[thing.sl[1]]))      
        
if savevtk is 'yes':
   t0=time.clock()
   results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
   t1=time.clock()
   t=t1-t0
   print 'time2writevtks',t
