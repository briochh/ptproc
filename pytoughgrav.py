# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 14:20:42 2014
\n
Module containing fuctions for building and running pytough models for gravity. \n
Defining grids\n
Restarting steady state models\n
Reading and plotting results - calculate gravity\n
Calculating relative gravity\n
Calculating and plotting rates of change\n
@author: Brioch Hemmings 
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import bisect
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import scipy.constants
import shutil
import time
import os

mpl.rcParams['xtick.labelsize']=14

mpl.rcParams['ytick.labelsize']=14

mpl.rcParams['axes.labelsize']=14
 
    
def grid1D( modelname, basefile, parameter, perm=5.0e-13, poro=0.34, rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}, cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]}, rechshift=0.0, wt=50.0 ):
    """ grid1D()\n
    Function to produce a 1D python grid.
    """
    print '\nbatch varying ' + parameter + '\n' \
    'creating ' + modelname + ' files with:\n' \
    'permeability = ' + str(perm) + ',\n' \
    'porosity = ' + str(poro) + ',\n' \
    'relative permeability dictionary =\n' \
    + repr(rp) + ',\n' \
    'capillary pressure dictionary =\n' \
    + repr(cp) + ',\n' \
    'recharge modifier = ' + str(rechshift) + '\n' \
    'and water table initial elevation ' + str(wt) + '.\n'   
    
    mod=modelname
    dat=basefile
    param=parameter
    
    
    dx=1
    dy=1
    #origin=[0,0,750]
    origin=[0,0,145]
    zcells=[1]*145 
    #print rechshift
    #zcells=[10]*34+[2]*80+[10]*19+[2]*30+[10]*25    
    
    geo = mulgrid().rectangular([dx],[dy],zcells, origin=origin, atmos_type =0, 
    convention = 2 )
    
    geo.atmosphere_volume= 1.0e50
    
    # write geometry to output file 
    geo.write(param+'/'+mod+'/2dgrd.dat') 
    
    ###### MAKE TOUGH GRID
    grid = t2grid().fromgeo(geo)
    
    # define relative permeability and cp paramters to use
    #
    #
    norp={'type':5, 'parameters':[]}
    #
    #
    nocp={'type':1, 'parameters':[0.0,0.0,1.0]}
    
    # define rock types and add cp and rp params
    lp=rocktype('lp   ', nad=3, permeability = [1.e-16]*2+[1e-16],
    porosity=0.1, conductivity=2.51, specific_heat=920) 
    lp.dry_conductivity=1.5 
    lp.tortuosity=0.0
    lp.relative_permeability=rp
    lp.capillarity=cp
    grid.add_rocktype(lp)
    
    hp=rocktype('hp   ', nad=3, permeability = [perm]*2+[perm],
    porosity=poro)
    hp.dry_conductivity=1.5 
    hp.tortuosity=0.0
    hp.relative_permeability=rp
    hp.capillarity=cp
    grid.add_rocktype(hp)
    
    b=rocktype('nocp ', nad=3, permeability = [perm]*2+[perm],
    porosity=poro)
    b.dry_conductivity=1.5 
    b.tortuosity=0.0
    b.relative_permeability=norp
    b.capillarity=nocp
    grid.add_rocktype(b)
    
    at=rocktype('atmos', nad=3, density=1.225, permeability = [perm]*2+[perm],
    porosity=1.0)
    at.dry_conductivity=1.5 
    at.tortuosity=0.0
    at.relative_permeability=norp
    at.capillarity=nocp
    grid.add_rocktype(at)
    
    
    # define rocktype of atmospher block    
    for blk in grid.atmosphere_blocks[:]:
        blk.rocktype= grid.rocktype['atmos']
        grid.block[(str(blk))].pmx=blk.rocktype.permeability[0]
        
    # assign rock properties
    # define low permeability region
    lam=0.004
#    k0=5.0e-13
    for blk in grid.blocklist[1:]:
        blk.rocktype = grid.rocktype['hp   ']
        # permeability modification
        col=geo.column[geo.column_name(str(blk))]
        lay=geo.column_surface_layer(col)    
        hmax=geo.block_surface(lay,col)
        pmx=blk.rocktype.permeability[0]*np.exp(-lam*(hmax-blk.centre[2]))
        grid.block[(str(blk))].pmx=pmx
        
    blay=geo.layerlist[-1]
    for col in geo.columnlist:
        blk=geo.block_name(blay.name,col.name)
        if blk in geo.block_name_list:
            grid.block[(blk)].volume=1E50
        
    # read template file    
       
    dat.parameter['print_block']='ee  1'
    # add rocktype, element and connection data to dat class
    dat.grid=grid    
    
    # INCON
    dat.incon.clear()
    # Define incon block
    initP=1.013e5
    initSG=0.99
    initT=25.0
    cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
    dat.incon[geo.block_name_list[0]]=cond
    for blk in grid.blocklist[1:]:
        if grid.block[np.str(blk)].rocktype==nocp:
            initP=1.013e5
            initSG=0.99
            initT=25.0
            cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
            dat.incon[np.str(blk)]=cond
        elif blk.centre[2] < wt:
            initP=1.013e5+(997.0479*9.81*np.abs(wt-blk.centre[2]))
            initSG=0.0
            initT=25.0
            cond=[[0.0,0.0,0.0],[initP,initSG,initT]]
            dat.incon[np.str(blk)]=cond
        elif grid.block[np.str(blk)].rocktype==lp:
            initP=1.013e5
            initSG=0.0
            initT=25.0
            cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
            dat.incon[np.str(blk)]=cond
        else:
            initP=1.013e5
            initSG=10.0
            initT=25.0
            cond=[[0.0,0.0,0.0],[1.013e5,initSG,initT]]
            dat.incon[np.str(blk)]=cond
           
    dat.generator.clear()
    
    # Define GENER block
#    fpms=7.7354e-6 # flux per meter squared
    fm=3.24e-8
    fc=-7.199e-7+((rechshift/10.)/3600/24)
    #print (rechshift/10.)/3600/24
    #print fc
    mingen=2.0e-7
    cols=[col for col in geo.columnlist]
    count=0
    dat.clear_generators()
    for col in cols:
        count=count+1
        lay=geo.column_surface_layer(col)
        blkname=geo.block_name(lay.name,col.name)
        gx=(grid.block[blkname].centre[2]*fm)+fc
        if gx < mingen: gx=mingen# for elevation dependant recharge!
        ex=1.0942e5
        gen=t2generator(name=' q'+col.name,block=blkname,type='COM1',gx=gx,ex=ex,hg=None,fg=None)
        #gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gx*col.area, ex=1.0942e5)
        dat.add_generator(gen) 
    
    dat.parameter['max_timestep']=3.0e6
    dat.parameter['print_interval']=30
    #dat.parameter['timestep']=[1000.0]
    #dat.output_times['time']=[1000.0,3600.0,8.6400e+04,3.1558e+07,3.1558e+08,3.1558e+09,3.1558e+10]
    #dat.output_times['num_times_specified']=7
    #dat.output_times['num_times']=7
    
    
           
    # write vtk of input information
    grid.write_vtk(geo,param+'/'+mod+'/inparam.vtk',wells=True)
       
    # write tough2 input file   
    dat.write(param+'/'+mod+'/flow2.inp')      
    shutil.copy('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/base/1D_20140620_2_py_it',param+'/'+mod+'/'+mod)


        
def readres( modelname, survey_points, save=False, savevtk=False, tough2_input=None, geom_data=None, results=None, fall=None):
    """ Function to read pytough results and calculate simulated changes in gravity associated with saturation changes.
    """
        ## read output file
    wells=survey_points
    mod=modelname

    
    t0=time.clock()
    if tough2_input is None:
       dat=t2data('flow2.inp') # tough2 input from input file
    else: dat=tough2_input

    if geom_data is None:
        geo=mulgrid('2dgrd.dat') # geometry from gempetry file
    else: geo=geom_data    

    if results is None:
        results=t2listing('flow2.out') # read output file
    
    grid=dat.grid # define input grid
    t1=time.clock()
    t=t1-t0
    print 'time2read .out=',t
    
    ## create directory for results
    if not os.path.exists('results'):
        os.makedirs('results')
    os.chdir('results')
    
    density=997.0479 # density of water
    yrsec=3600*24*365.25 # 1 year in seconds
    
    
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
    #ntheta=1000.0 # number of increments
    #dtheta=2.0*np.pi/ntheta # length of element of integral
    #thetas=np.arange(0.0,2.0*np.pi,dtheta) #array of integral element lengths
    
    
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
        if geo.atmosphere_type is 2: 
            stpoint=0
        else: stpoint=1
        for lay in geo.layerlist[stpoint:]:
            #print "lay.name =", lay.name
            #print "col.name =",col.name
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
        while results.time/yrsec <= 110 and count < results.times.size:
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
        microgal=well_water_mass*2*np.pi*scipy.constants.G*(10**8)/col.area   
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
        intgravplt=plt.figure()
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
        
        if save:  
           t0=time.clock() 
           if not os.path.exists('mark2'):
               os.makedirs('mark2')
               os.chdir('mark2')
           zt_density_matrix=np.concatenate((
           [np.concatenate((np.array([0]),times))],
            np.concatenate((zlist.reshape(len(zlist),1),wellro),
                           axis=1)),
                           axis=0)
           if fall is not None: 
               fall.write(str(mod)+'\t'+str(well_water_mass.min())+'\t'+str(well_water_mass.max())+'\t'+str(well_water_mass.max()-well_water_mass.min())+'\t'+str(microgal.min())+'\t'+str(microgal.max())+'\t'+str(microgal.max()-microgal.min())+'\n')
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
           np.savetxt('ztro.dat',zt_density_matrix)                
           np.savetxt('waterweight'+str(wellno)+'.dat',zip(times,well_water_mass))
           np.savetxt('microgal'+str(wellno)+'.dat',zip(times,microgal))
           np.savetxt('axsym_int_microgal'+str(wellno)+'.dat',zip(times,(gravt-gravt[0])))
           im.savefig(mod+'_elemental_cont'+str(wellno)+'.pdf')       
           im2.savefig(mod+'_sl_t_profile'+str(wellno)+'test.png',dpi=300)
           #im2.savefig('sl_t_profile'+str(wellno)+'.eps')
           im1.savefig(mod+'_boug'+str(wellno)+'.pdf')
           intgravplt.savefig(mod+'_axsym_int_grav'+str(wellno)+'.pdf')
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
            
    if savevtk:
       t0=time.clock()
       results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
       t1=time.clock()
       t=t1-t0
       print 'time2writevtks',t

def grate( modelname, in_ts, winlen=[2,5,10], save=True, input_in="yrs", fall=None, fallmax=None ):
    """ grate( timeseries, window_length, save_option, input_in )\n
    Calculate and plot simulated rates of gravity change. 
    """
    
    mod=modelname
    print 'window lengths =',winlen
    if save: 
        print 'save is on'
    else:
        print 'save is off'
    print 'input timeseries assumed to be in',input_in
    
    plottimes={}
    windg={}
    yrsec=3600*24*365.25
    if input_in is 'yrs': in_ts[:,0]=in_ts[:,0]*yrsec
    for win in winlen:
        n=0
        plott=np.empty(len(in_ts))
        dy=np.empty(len(in_ts))
        yrdg=np.empty(len(in_ts)) 
        for t,g in in_ts:
            te=t+win*yrsec 
            i=bisect.bisect(in_ts[:,0], te)
            if i == len(in_ts):
                dy[n] = in_ts[i-1,1]+((te-in_ts[i-1,0])*(in_ts[i-1,1]-in_ts[i-2,1])/(in_ts[i-1,0]-in_ts[i-2,0]))
            else:
                dy[n]=in_ts[i-1,1]+((in_ts[i,1]-in_ts[i-1,1])*(te-in_ts[i-1,0])/(in_ts[i,0]-in_ts[i-1,0]))
            yrdg[n]=dy[n]-g
            plott[n]=(t+win*yrsec/2)
            n+=1
        plottimes['win_'+str(win)]=plott
        windg['win_'+str(win)]=yrdg
    
    #times=plott/yrsec
        
    #im1=plt.figure()
    #plt.plot(times,yrdg)
    #plt.ylabel(r'$\Delta g$ (microgal)/'+str(win)+'yr')
    #plt.xlabel('Time (years)')
    #plt.axis([0.0, 110,None,None])
    #
    #im1.savefig('mugal_per_'+str(win)+'yr.pdf')  
    
    x = in_ts[:,0]/yrsec
    y = in_ts[:,1]
    f = interp1d(x, y)
    #f2 = interp1d(x, y, kind='cubic')
    
    
    
    xnew = np.linspace(x.min(), x.max(), 100000)
    ynew=f(xnew)
    dgbydt=np.gradient(ynew,np.diff(xnew)[0])
    
    output=np.array([['peryear',dgbydt.min(),dgbydt.max()]])
    
    im2, axarr=plt.subplots(2,sharex=True)
    axarr[0].tick_params(axis='both',labelsize=18)
    axarr[1].tick_params(axis='both',labelsize=18)
    
    data, =axarr[0].plot(x,y,'-',label='data',linewidth=2)
    #interp, = plt.plot(xnew,f(xnew),'b-',label='linear')
    axarr[0].set_ylabel(r'$\Delta g$ ($\mu$gal)',fontsize=18)
    #plt.xlabel('Time (years)') 
    #ax2=plt.twinx()
    color_cycle=axarr[1]._get_lines.color_cycle
    next(color_cycle)
    #grad, = axarr[1].plot(xnew,dgbydt,'-',label='gradient', color='0.6' )
    #plt.legend(loc='best')
    leghand=[data]#,grad]
    leglab=['data']#,'$\Delta g$/yr']
    for win in winlen:
        temp=axarr[1].plot( plottimes['win_'+str(win)]/yrsec,windg['win_'+str(win)],'-', markersize=12,label=r'$\Delta g$/'+str(win)+'yrs',linewidth=2)
        leghand=leghand+temp
        leglab=leglab+[r'$\Delta g$/'+str(win)+'yrs']
        output=np.concatenate((output,[[win,windg['win_'+str(win)].min(),windg['win_'+str(win)].max()]]))
    axarr[1].set_ylabel(r'$\Delta g$ per time ($\mu$gal/time)',fontsize=18)
    plt.xlabel('Time (years)',fontsize=18)    
    data.axes.legend(leghand,leglab,bbox_to_anchor=(0., -0.05, 1., .202), loc='upper center',
           ncol=4, mode="expand", borderaxespad=0.,fontsize=18)
    plt.axis([0.0, 100,None,None])
    
    plt.show()
    
    
    if save:
           im2.savefig(mod+'_mugal_all_per_yr.pdf')
           np.savetxt(mod+'_minmaxs.txt',output,fmt='%s')
           if fall is not None:           
               fall.write('modelname = '+str(mod)+ '\n'
               'byyear \t gradmin \t gradmax \n')
               np.savetxt(fall,output,fmt='%s',delimiter='\t', newline='\n')
               fall.write('\n')
           if fallmax is not None:
               fallmax.write(str(mod)+'\t'+str(np.max(np.abs([float(i) for i in output[0][1:]])))+'\t'+str(np.max(np.abs([float(i) for i in output[1][1:]])))+'\t'+str(np.max(np.abs([float(i) for i in output[2][1:]])))+'\t'+str(np.max(np.abs([float(i) for i in output[3][1:]])))+'\n')
           #f = open('resultxt_'+str(wellno)+'test.txt','w')
    return output               

def relgrav(reference_modelname,reference_ts,test_modelname,test_ts):
    ref_mod=reference_modelname
    mod=test_modelname  
    ref_grav=reference_ts
    grav=test_ts

