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
#import matplotlib.mlab as ml
import numpy as np
import bisect
from scipy.interpolate import interp1d
#from scipy.interpolate import griddata
import scipy.integrate as integrate
import scipy.constants
import shutil
import time
import os
import copy
import numpy.ma as ma
import cPickle as pickle
import random

mpl.rcParams['xtick.labelsize']=14

mpl.rcParams['ytick.labelsize']=14

mpl.rcParams['axes.labelsize']=14
 
def save_obj(obj, name ):
    with open( name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name , 'rb') as f:
        return pickle.load(f)
    
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

def geo2D( modelname, length = 500., depth = 500., width = 1., celldim = 10.,
          origin = ([0,0,0]), xcells = None,  zcells = None,
          surface = None, wells = None, atmos_type=0, min_thick=2.0):
    """ Function to generate a 2D grid tough2 grid."""
    mod=modelname
    if xcells is None:
        xcells=[celldim]*int((length/celldim))
        
    if zcells is None:
        zcells=[celldim]*int((depth/celldim))
    dy=[width]
    #zcells=[10]*34+[2]*80+[10]*19+[2]*30+[10]*25
    geo = mulgrid().rectangular(xcells,dy,zcells, origin=origin, atmos_type =atmos_type, 
    convention = 2 )  # creates geometry 20 cells that are 500 m width in x,
    # 1 cell 1000 m width in y
    # 20 cells 100 m high in z  
    # to make use of more possible grid names add char=ascii_lowercase+ascii_uppercase
    geo.atmosphere_volume= 1.0e50 # change volume of atmos cell to 1e50
    if wells is not None:
        colstorefine=[geo.columns_in_polygon([np.concatenate((np.add(well,(-15,0)),[origin[2]-sum(zcells)])),np.concatenate((np.add(well,(15,0)),[origin[2]]))]) for well in wells]
        geo.refine(np.hstack(colstorefine),bisect='y',chars = '123456789')
        colstorefine=[geo.columns_in_polygon([np.concatenate((np.add(well,(-5,0)),[origin[2]-sum(zcells)])),np.concatenate((np.add(well,(5,0)),[origin[2]]))]) for well in wells]
        geo.refine(np.hstack(colstorefine),bisect='y',chars = '123456789')
        colstorefine=[geo.column_containing_point(well) for well in wells]
        geo.refine(np.hstack(colstorefine),bisect='y',chars = '123456789')
    if surface is not None:
        geo.fit_surface(surface, silent=True, layer_snap=min_thick) # fit topograpghy surface
    geo.write(mod+'/grd.dat')
    # write geometry to output file 
    return geo
    
def grid2D(modelname,geo,dat,rocks,boundcol,lpregion=[[0,0,0],[0,0,0]],satelev=0.0,atmosP=1.013e5,pmx_lamda=0.004):
    """ Function to add TOUGH2 data to geometry"""
    grid = t2grid().fromgeo(geo)
    atmos=grid.atmosphere_blocks
    for rock in rocks:
        grid.add_rocktype(rock) 
    # loop over every element block in the grid
    for blk in grid.blocklist[0:]:
        if blk not in atmos:
            col=geo.column[geo.column_name(str(blk))] # column containing current block
            tlay=geo.column_surface_layer(col)    
            hmax=geo.block_surface(tlay,col)
            lay=geo.layer[geo.layer_name(str(blk))] # layer containing current block
            if col in boundcol: # if block is in lateral boundary columns
                rocktype='nocp '
                pmx=grid.rocktype[rocktype].permeability[0]*np.exp(-pmx_lamda*(hmax-blk.centre[2])) # calculating depth dependent permeability modifier
                initP=atmosP
                initSG=0.99
                initT=25.0
                #send to function to assign a blocktype and initial condition and pmx.
                rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx,infvol=True)
            elif (blk.centre[2] > lpregion[0][2] and 
                blk.centre[2] <= lpregion[1][2] and 
                blk.centre[0] > lpregion[0][0] and 
                blk.centre[0] <= lpregion[1][0]):
                rocktype='lp   '
                pmx=grid.rocktype[rocktype].permeability[0]
                initP=atmosP
                initSG=10.7
                initT=25.0                
                rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx)                  
            else:
                rocktype='hp   '
                pmx=grid.rocktype[rocktype].permeability[0]*np.exp(-pmx_lamda*(hmax-blk.centre[2]))
                initP=atmosP
                initSG=10.7
                initT=25.0                
                rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx)
            if blk.centre[2] < satelev:
                initP=1.013e5+(997.0479*9.81*abs(blk.centre[2]))
                initSG=0.0
                initT=25.0
                pmx=None # dont change pmx
                rockandincon(blk,grid,dat,None,initP,initSG,initT,pmx)
        else:
            rocktype='atmos'
            pmx=grid.rocktype[rocktype].permeability[0]          
            initP=atmosP
            initSG=0.99 # initial gas saturation  
            initT=25.0 # initial temperature - TOUGH2 doesn't seem to like < 1.0 C
            rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx)
    return grid
        
def rockandincon(blk,grid,dat,rocktype,P,SG,T,pmx,infvol=False):   
    if rocktype is not None:
        blk.rocktype=grid.rocktype[rocktype]
    dat.incon[str(blk)]=[None,[P,SG,T]]
    if pmx is not None:
        grid.block[(str(blk))].pmx=pmx
    if infvol:
        grid.block[str(blk)].volume=grid.block[str(blk)].volume*1E50 # more robust later on as we retain somthing of orginal volume
#       grid.block[str(blk)].volume=1E50 
        
def topsurf(surfpath,delim='\t',headerlines=1,width=10):
    """ reads and reshapes surface profile for use in 2D model """
    ## top surface
    surf = np.loadtxt(surfpath,delimiter=delim,skiprows=headerlines)
    #np.loadtxt(r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\2Ddev\2dprof.txt', delimiter='\t', skiprows=1) # load surface file
    halfwidth=width/2.0
    neghalf=-halfwidth
    
    minw=neghalf*np.ones((surf.shape[0],1)) # adapt to min max of y (-5,+5)
    maxw=halfwidth*np.ones((surf.shape[0],1))
    surf=np.concatenate(((np.concatenate((
    np.hsplit(surf,2)[0],minw,np.hsplit(surf,2)[1]),axis=1)),
    (np.concatenate((np.hsplit(surf,2)[0],maxw,np.hsplit(surf,2)[1]),axis=1))),
    axis=0)
    return surf
def makeradial(geo,grid,width=1.):
    """turn 2D grid into radial grid about x=0"""
    if grid is not None:
        if len(grid.atmosphere_blocks) == 1:
            stpoint=1
        else:
            stpoint=0
        for blk in grid.blocklist[stpoint:]:
            blk.volume=blk.volume*2.*np.pi*blk.centre[0]/width
            #col=geo.column[geo.column_name(str(blk))] # column containing current block
        for conn in grid.connection.values():
            #print conn
            #print conn.direction
            if conn.direction==3:
                R=conn.block[0].centre[0]
            elif conn.direction==1:
                cellRd=zip([blk.centre[0] for blk in conn.block],[cdist for cdist in conn.distance])
                cellRd.sort()
                R=sum(cellRd[0])
            conn.area=2*np.pi*R*conn.area/width
    
    for col in geo.columnlist:
        col.area=2*np.pi*col.centre[0]*col.area/width
    geo.radial=True
           
def gen_constant(mod,geo,grid,dat,constant=7.7354e-6,elev_m=None,elev_c=None,mingen=2.0e-7,enthalpy=1.0942e5):
    dat.clear_generators()
    f = open(mod+'/genertot.txt','w')
    f.write('Model = '+mod+'\n')
    allgens=[]
    cols=[col for col in geo.columnlist]
    if elev_m is None:
        f.write('Constant generation ='+str(constant)+' kg/s/m2\n')
        for col in cols:
            lay=geo.column_surface_layer(col)
            blkname=geo.block_name(lay.name,col.name)
            gx=constant
            gxa=col.area*gx
            gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gxa, ex=enthalpy)
            dat.add_generator(gen)
            allgens.append(gxa)
    else:
        f.write('Elevation dependent generation \n'
                'gen =' +str(elev_m)+ '*z +' +str(elev_c)+ '\n')
        for col in cols:
            lay=geo.column_surface_layer(col)
            blkname=geo.block_name(lay.name,col.name)
            gx=(grid.block[blkname].centre[2]*elev_m)+elev_c
            if gx < mingen:
                gx=mingen
            gxa=col.area*gx
            gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gxa, ex=enthalpy)
            dat.add_generator(gen)
            allgens.append(gxa)
    allgens=np.array(allgens)
    gensum=np.sum(allgens)
    f.write('Total generation in model = '+str(gensum)+' kg/s\n')
    f.write('Total generation rate per m2 = '+str(gensum/geo.area)+' kg/s\n')
    f.close()
    
def gen_variable(mod,geo,grid,dat,ts="C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/2Ddev/rand.dat",season_bias=0.65,length=100,
                 wavelength=1,maxlength=3e5,new_rand=None,constant=7.7354e-6,
                 elev_m=None,elev_c=None,mingen=2.0e-7,enthalpy=1.0942e5):
    """define time dependent generation rate for recharge"""
    dat.clear_generators()
    yrsec=3600*24*365.25
    wavelength=wavelength*yrsec
    length=length*yrsec
    maxlength=maxlength*yrsec
    fm=elev_m
    fc=elev_c
    mult=season_bias
    knownts=False    
    
    allgens=[]
    if new_rand<>None:
        ts=[]
        times=[0.0]+np.arange((wavelength/2),length,wavelength/2).tolist()+[length,maxlength]
        numt=len(times)
        for i in xrange(0,numt-3,2): ts=ts+[random.gauss(1,new_rand)]
        np.savetxt(mod+'/rand.dat',ts)
    elif isinstance(ts,basestring) and os.path.isfile(ts): 
        ts=np.loadtxt(ts) # load random data file
        np.savetxt(mod+'/rand.dat',ts)
    elif type(ts).__module__ ==  np.__name__:
        if np.shape(ts.shape) == 2: # if timeseries is complete with time info
            knownts=True
            np.savetxt(mod+'/init_rech.dat',ts)
    #print type(ts)
    if knownts:
        #print knownts
        times=[0.0]+np.cumsum(ts[:,1]).tolist()+[maxlength]
    else:
        times=[0.0]+np.arange((wavelength/2),length,wavelength/2).tolist()+[length,maxlength]
    #print times
    numt=len(times)
    xs=[]    
    zs=[]
    Areas=[]
    for col in geo.columnlist:
        gxc=[]
        lay=geo.column_surface_layer(col)
        blkname=geo.block_name(lay.name,col.name)
        xs=xs+[grid.block[blkname].centre[0]]
        zs=zs+[grid.block[blkname].centre[2]]
        if elev_m is None:
            gx=constant
        else:
            gx=(grid.block[blkname].centre[2]*fm)+fc
        if gx < mingen: gx=mingen
        # for elevation dependant recharge!
        if knownts:
            for month in ts:
                gxc=gxc+[((grid.block[blkname].centre[2]*fm)+(fc))*month[0]]
        else:
            for i in xrange(0,numt-3,2):
                highgx=((1+mult)*((grid.block[blkname].centre[2]*fm)+(fc)))*ts[i/2]
                if highgx < (1+mult)*mingen: highgx=(1+mult)*mingen
                lowgx=((1-mult)*((grid.block[blkname].centre[2]*fm)+(fc)))*ts[i/2]
                if lowgx < (1-mult)*mingen: lowgx=(1-mult)*mingen
                gxc=gxc+[lowgx,highgx]
    
        gxc=gxc+[gx,gx]
        ex=numt*[1.0942e5]
        Areas=Areas+[col.area]
        gxa=np.multiply(col.area,gxc).tolist()
        allgens.append(gxa)
        gen=t2generator(name=' q'+col.name,block=blkname,type='COM1',gx=None,ex=None,hg=None,fg=None, rate=gxa, enthalpy=ex,   time=times,ltab=numt,itab=numt-1)
        #gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gx*col.area, ex=1.0942e5)
        dat.add_generator(gen)
    
    dat.output_times['time_increment']=2.4192E6
    dat.output_times['time']=[1.0]+times[1:100]
    dat.output_times['num_times_specified']=len(dat.output_times['time'])
    dat.output_times['num_times']=200
    allgens=np.array(allgens)
    gensum=np.sum(allgens,axis=0) # <<< new >>>  old >>> sum(row[:] for row in allgens)
    tforplot=[times[0]]
    tforplot=np.append(tforplot,[2*[j] for j in times[1:-1]]+[times[-2]+yrsec])
    tforplot=np.hstack(tforplot)
    gforplot=[2*[j] for j in gensum[0:-1]]
    gforplot=np.hstack(gforplot)
    Area=sum(Areas)
    fig,ax1=plt.subplots()
    ax1.plot(tforplot/yrsec,gforplot/Area)
    ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
    ax1.set_ylabel(r'Generation rate (kg/s/m$^2$)')
    ax1.set_xlabel('Time (years)')
    ax2=ax1.twinx()
    print max(gforplot/Area)
    ax2.plot()
    ax2.set_ylim(0,max(gforplot/Area)*3600*24)   
    ax2.set_ylabel(r'Equivalent recharge rate (mm/d)')
    
    plt.savefig(mod+'/rech.pdf')
    np.savetxt(mod+'/genertot.txt',np.vstack((tforplot,gforplot)).T)
    return allgens,xs,zs,Areas,times
        
def readres( modelname, survey_points, save=False, savevtk=False, tough2_input=None, geom_data=None, results=None, sat={}, fall=None, maxtime=20):
    """ Function to read pytough results and calculate simulated changes in gravity associated with saturation changes.
    """
    ## read output file
    wells=survey_points
    mod=modelname

    
    t0=time.clock()
    if type(tough2_input) is not t2data and tough2_input is not None:
        raise TypeError('data needs to be type t2data. Type ' + str(type(tough2_input)) + ' found. You idiot')
    elif tough2_input is None:
        raise TypeError('data needs to be type t2data. Currently None found. You idiot')
    else: dat=tough2_input
    if type(geom_data) is not mulgrid and geom_data is not None:
        raise TypeError('data needs to be type mulgrid. Type ' + str(type(geom_data)) + ' found. You idiot')
    elif geom_data is None:
        raise TypeError('data needs to be type mulgrid. Currently None found. You idiot')    
    else: geo=geom_data
    if type(results) is not t2listing and results is not None:
        raise TypeError('results needs to be type t2listing. Type ' + str(type(results)) + ' found. You idiot')
    elif results is None:
        print('No Results files (flow.out) passed. please read flow2.out and pass to ptg.readres')
        

    
    grid=dat.grid # define input grid
    width=geo.bounds[1][1]-geo.bounds[0][1]   
    if not geo.radial:    
        makeradial(geo,None,width=width)

    
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
        if blk.volume >= 1.0E50:
            dz=np.concatenate((dz,[blk.volume/geo.column[geo.column_name(blk.name)].area/1.0E50]))
            #dz=np.concatenate((dz,[blk.volume/geo.column[geo.column_name(blk.name)].area/1.0E47]))
        else: dz=np.concatenate((dz,[blk.volume/geo.column[geo.column_name(blk.name)].area])) # array of thickness of each block
        # dz=np.concatenate((dz,[geo.layer[geo.layer_name(blk.name)].thickness])) # array of thickness of each block
        dx=np.concatenate((dx,[geo.column[geo.column_name(blk.name)].side_lengths[0]])) # array of x direction cell widths
        cen=np.concatenate((cen,[blk.centre]))
    
    xzarea=dz*dx # array of cell x by z area
    xs=np.array([c[0] for c in cen])
    zs=np.array([c[2] for c in cen])
    xli=np.unique(xs)
    zli=np.unique(zs)
    xi,zi=np.meshgrid(xli,zli)
    index=[]
    jndex=[]
    for xz in zip(xs,zs):
        index.append(np.where(xz[0]==xi[0])[0][0])
        jndex.append(np.where(xz[1]==zi[:,0])[0][0])
        
        
    ## numerical solution of ring integral from 0 to 2pi slightly quicker
    ## set up arrays of theta increments
    #ntheta=1000.0 # number of increments
    #dtheta=2.0*np.pi/ntheta # length of element of integral
    #thetas=np.arange(0.0,2.0*np.pi,dtheta) #array of integral element lengths
    t1=time.clock()    
    print 'time to initialise results grid=',(t1-t0)
    
    #%%
    ## loop over each gravity survey point given in wellx
    wellno=0
    #maxtime=maxtime
    outtimes=[outtime for outtime in results.times if outtime/yrsec <= maxtime]
    #sat={}  # set or resest dictionary of saturations
    
    for well in wells:
        twell=time.clock()
        wellno=wellno+1
        print well
        wellblk=[] # list for block names in well column
        zlist=[] # array of cell depths in column
        
        
        ## generate list of cells beneath survey point (within well)
        ## usedfor bouguer slab
        col=geo.column_containing_point(well) # column that contains the survey point
        if geo.atmosphere_type is 2: 
            stpoint=0
        else: stpoint=1
        for lay in geo.layerlist[stpoint:]:
            if geo.block_name(lay.name,col.name) in geo.block_name_list:
                blk=grid.block[geo.block_name(lay.name,col.name)]
                wellblk.append(blk) # list of cells in well column (used in bouguer slab calculation)
                zlist=np.concatenate((zlist,[blk.centre[2]]))
        
        t1=time.clock()
        t=t1-twell
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
            # distance to every element from survey point - Theta must be set above
            # s=((alpha**2)-(2*alpha*well[0]*cos(thetas))+(well[0]**2)+((zp-blkz)**2))**(1./2.) 
            # val=(sum(dtheta/s**3)) # fractional contribution of current ring
            # python integrate method: (slower)?
            I= lambda theta: 1/(((alpha**2)-(2*alpha*well[0]*np.cos(theta))+(well[0]**2)+((zp-blkz)**2))**(1.5)) 
            val,err= integrate.quad(I, 0.0, 2*np.pi)
            # multiply by radius of ring, vertical distance from survey point to ring, xzarea 
            comp.append((zp-blkz)*alpha*val*blkxzarea)
            # comp.append((zp-blkz)*val*blkxzarea)
            i+=1
        t2=time.clock()
        print 'time to calculate contribution',(t2-t1)
            
        
        # loop over results table untill desired time (in years)
        count=0
        for outtime in outtimes:
            if count < results.times.size:
                t0=time.clock()
                print('On Station %d of %d' % (wellno,len(wells)))
                print('timestep %d out of %d (%5.2f yrs)' % (count+1,results.times.size,outtime/yrsec))
                if str(outtime) not in sat: # for first well pull out saturation for each tstep                
                    sat[str(outtime)]=copy.copy(results.element['SL'][1:]) # pull out saturation index [0] is the atmosphere cell so no use to us                
                    results.next() # move to next timestep             
                dg=[] # array for collecting together elemental contributions of integral method
                twellsl=[] # array for collecting saturation in well for current timestep
                twellrho=[] # array for collecting saturation density in well for current timestep
                twell_water_mass=0 # starting point for water mass in well at current timestep
                # loop over every block for current timestep
                i=0
                for blk in grid.blocklist[1:]: # dont bother with atmos cell
                    blkrho=sat[str(outtime)][i]*blk.rocktype.porosity*density # saturation density in element
                    dg.append(comp[i]*blkrho) # gravity conribution of element 
            
                    ## For bouguer slab.....
                    #calculated saturation, saturation density and water mass in column at current timestep
                    if blk in wellblk:
                        #outd[blk.name].append(sat[i])
                        twellsl.append([sat[str(outtime)][i]])
                        twellrho.append([blkrho])
                        #if blk is not wellblk[-1]:
                        #wellvol=wellvol+blk.volume
                        col=geo.column[geo.column_name(str(blk))]
                        if blk.volume >= 1.0E50:
                            dummyvolume=width*np.abs(col.bounding_box[1][0]-col.bounding_box[0][0])*((blk.volume/1.0E50)/col.area)
                        else: dummyvolume=width*np.abs(col.bounding_box[1][0]-col.bounding_box[0][0])*(blk.volume/col.area)                        
                        #print geo.column[geo.column_name(str(blk))].area
                        #print 'dummyvolume=',dummyvolume
                        twell_water_mass=twell_water_mass+(blkrho*dummyvolume)
                        # wellsl=np.concatenate((wellsl,[outsl1[1]]))
                        # wellro=np.concatenate((wellro,[outsl1[1]*density*blk.rocktype.porosity]))
                    i+=1 # inrement to next element
                #compile array of saturation, density and water mass in current column for all timesteps
                wellsl=np.concatenate((wellsl,twellsl),axis=1)
                wellro=np.concatenate((wellro,twellrho),axis=1)
                well_water_mass.append(twell_water_mass)
                
                # Total axissummetric intgral gravity due to water mass in model for current time step
                grav=6.67e-3*sum(dg)
                gravt.append(grav) # compile gravity due to water at current survey point for each timestep 
                print grav
                count+=1 # incrment timestep counter
                dt=time.clock()-t0
                print('time for timestep = %e s' % dt)

        # Bouguer slab gravity approximations
        well_water_mass=np.array(well_water_mass)
        #print ((col.area*width)/(2*np.pi*col.centre[0]))
        microgal=well_water_mass*2*np.pi*scipy.constants.G*(10**8)/((col.area*width)/(2*np.pi*col.centre[0])) # bit of trickery to get 2D distribution....   
        microgal=microgal-microgal[0] # gravity difference    
        
        t1=time.clock()
        t=t1-twell
        print 'time4well_calculations',t    
        #%%
        ## Plot results for current well
        t0=time.clock()
        times=results.times[0:count]/yrsec # convert times calculted to yrs 
        #gcont=ml.griddata(xs,zs,np.array(dg/xzarea)*6.67e-3,xi,zi,interp='linear')
        c=0
        dumgrid=np.empty(xi.shape)*np.NaN
        dumgrid2=np.empty(xi.shape)*np.NaN
        for j,i,g,a in zip(jndex,index,np.array(dg)*6.67e-3,xzarea):
            dumgrid[j,i]=g/a
            dumgrid2[j,i]=100.*(g/grav)
            c=c+1
        gcont=ma.array(dumgrid,mask=np.isnan(dumgrid))        
        gcont2=ma.array(dumgrid2,mask=np.isnan(dumgrid2))   
        #gcont=griddata(np.vstack((xs,zs)).T,np.array(dg/xzarea)*6.67e-3,(xi,zi),method='nearest')
        elcont=plt.figure(figsize=[8,3.6]) 
        plt.pcolormesh(xi,zi,gcont,shading='flat',edgecolor='face', rasterized=True)
        plt.colorbar().set_label(r'Contribution to g (' + r'$\mu$'+'gal/m' + r'$^{2}$'+')')
        plt.xlim((xs.min(),xs.max()))
        plt.ylim((zs.min(),zs.max()))
        
        elcont2=plt.figure(figsize=[8,3.6]) 
        plt.pcolormesh(xi,zi,gcont2,shading='flat',edgecolor='face', rasterized=True)
        plt.colorbar().set_label('Percent contribution of element\n to total water induced gravity (%)')
        plt.xlim((xs.min(),xs.max()))
        plt.ylim((zs.min(),zs.max()))
        # test plot of contibutions
#        im=plt.figure(figsize=[8,3.6])
#        plt.scatter(xs, zs, c=np.array(dg/xzarea)*6.67e-3, edgecolor='none', marker='s')
#        plt.colorbar()
#        plt.xlim((xs.min(),xs.max()))
#        plt.ylim((zs.min(),zs.max()))
#        plt.show()
        
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
        profplt=plt.pcolormesh(T,Z,wellsl,cmap=cm.jet_r,vmin=0.0,vmax=0.8,shading='flat',rasterized=True)
        plt.axis([T.min(), T.max(), 0, Z.max()])
        plt.ylabel('Z (m)')
        plt.xlabel('Time (years)')    
        cbar=plt.colorbar(profplt,orientation='vertical')
        cbar.set_label('Saturation')   
        t1=time.clock()
        t=t1-t0
        print 'time2plotwell',t
        #%%
        if save:  
           t0=time.clock() 
           #if not os.path.exists('mark2'):
           #    os.makedirs('mark2')
           #    os.chdir('mark2')
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
           elcont.savefig(mod+'_elemental_cont'+str(wellno)+'.pdf',dpi=400) 
           elcont2.savefig(mod+'_elemental_cont_norm'+str(wellno)+'.pdf',dpi=400)
           im2.savefig(mod+'_sl_t_profile'+str(wellno)+'.pdf',dpi=400)
           #im2.savefig('sl_t_profile'+str(wellno)+'.eps')
           im1.savefig(mod+'_boug'+str(wellno)+'.pdf')
           intgravplt.savefig(mod+'_axsym_int_grav'+str(wellno)+'.pdf')
           #
           #savefig('sl_t_profile.pdf')
           t1=time.clock()
           t=t1-t0
           print 'time2saveplot',t
           plt.close('all')
        
        t1=time.clock()
        t=t1-twell
        print 'total time for well',t
        
            
    if savevtk:
       t0=time.clock()
       results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
       t1=time.clock()
       t=t1-t0
       print 'time2writevtks',t
    if save:
        if os.path.isfile('sat.pkl'):
            print('sat alreay pickled')
        else:
            save_obj(sat,'sat.pkl')
    return results,sat
#%%
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

def relgrav(reference_modelname,test_modelname,reference_ts='axsym_int_microgal5.dat',test_ts=['axsym_int_microgal1.dat'],save=True,time_in='yrs'):
    plt.close('all')    
    cd=os.getcwd()   
    ref_mod=reference_modelname
    mod=test_modelname  
    #ref_grav=reference_ts
    
    os.chdir(ref_mod+'/results')
    ref_grav=np.vstack(np.loadtxt(reference_ts)).T
    
    os.chdir(cd)
    if type(mod) is str:
        if type(test_ts) is not str:
            modlist=[mod]*len(test_ts)
        else:
            test_ts=[test_ts]
            modlist=[mod]
            
    num=1        
    for ts,mod in zip(test_ts,modlist):
        os.chdir(cd)
        os.chdir(mod+'/results')
        grav=np.vstack(np.loadtxt(ts)).T
        if grav[0][-1] > ref_grav[0][-1]:
            ref_grav=np.concatenate((ref_grav,np.array([[grav[0][-1]],[grav[1][-1]]])),axis=1)
        im1=plt.figure()
        if time_in is 'yrs':
            yrsec=1
        else: yrsec=3600*24*365.25
        plt.plot(ref_grav[0]/yrsec,ref_grav[1],'-',linewidth=2,label='Reference Station')
        plt.plot(grav[0]/yrsec,grav[1],'-',linewidth=2,label='Benchmark P'+str(num))


        times=np.array([])
        gravdif=np.array([])
        i=0
        for t,g in zip(grav[0],grav[1]):
           #print 't='+str(t)
           #print 'g='+str(g)
           while i < len(ref_grav[0]):
              #print i
              if t >=ref_grav[0][i] and t <= ref_grav[0][i+1]:
                 times=np.concatenate((times,[t]),axis=0)
                 #print 'times='+str(times)
                 grad=((ref_grav[1][i+1])-(ref_grav[1][i]))/((ref_grav[0][i+1])-(ref_grav[0][i]))
                 # print 'grad='+str(grad)
                 refgatt=ref_grav[1][i]+((t-ref_grav[0][i])*grad)
                 #print 'refgatt='+str(refgatt)
                 gravdif=np.concatenate((gravdif,[g-refgatt])
                 ,axis=0)
                 break
              else: i=i+1    
                    
                
        plt.plot(times/yrsec,gravdif,linewidth=2,label='Relative gravity') 
        plt.xlim((0,times.max()/yrsec))
        plt.ylim((-70,70))
        plt.ylabel(r'$\Delta g$ (microgal)',fontsize=18)
        plt.xlabel('Time (years)',fontsize=18)  
        plt.tick_params(axis='both',labelsize=18)
        plt.legend()
        if save:                 
           np.savetxt('gravdiff'+str(num)+'.dat',zip(times,gravdif))
           im1.savefig('vsref'+str(num)+'.pdf')   
           im1.savefig('vsref'+str(num)+'.png')
           f = open('relresultxt'+str(num)+'.txt','w')
           if mod is not ref_mod:          
               f.write('Model = '+mod+'\n'
                           'Base = '+ref_mod+'\n'
                           'relgrav max (microgal) =' +str(gravdif.max())+'\n'
                           'relgrav min (microgal) =' +str(gravdif.min())+'\n'
                           'Max amplidute (relgrav)='+str(gravdif.max()-gravdif.min())+'\n')
           else:
               f.write('Model = '+mod+ts+'\n'
                       'Base = '+ref_mod+reference_ts+'\n'
                       'relgrav max (microgal) =' +str(gravdif.max())+'\n'
                       'relgrav min (microgal) =' +str(gravdif.min())+'\n'
                       'Max amplidute (relgrav)='+str(gravdif.max()-gravdif.min())+'\n')
               
           f.close()
        num=num+1
    os.chdir(cd)

