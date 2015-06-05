# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 09:57:38 2015
\n
Module containing fuctions for building and running pytough icey models. \n

Reading results\n

@author: Brioch Hemmings 
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import numpy as np
import scipy as spy
import shutil
import time
import os
import pytoughgrav as ptg

def icegeo( modelname, length = 500., depth = 500., width = 1., celldim = 10.,
          origin = ([0,0,0]), xcells = None,  zcells = None,
          surface = None, wells = None, atmos_type=0, min_thick=2.0):
    """ 
    method for defining geometry of ice work. 
    
    Initially it is the same as ptg definition so is just passed to that.
    
    """
    geo=ptg.geo2D( modelname, width=width, celldim=celldim, origin=origin,
              zcells=zcells, xcells=xcells, surface=surface, atmos_type=atmos_type, min_thick=min_thick)
    return geo

def icegrid(geo,dat,rocks,boundcol,lpregion=None,hpregion=None,heatsource=None,satelev=0.0,atmosP=1.013e5,pmx_lamda=0.004, glacier_limit=2500.):
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
            Tmin=2      
            if blk.centre[0] <= glacier_limit: 
                initT=Tmin # initial temperature - TOUGH2 doesn't seem to like < 1.0 C
            else:
                initT = 25.8 - (hmax*(5.4/1000)) # 15.+((2000.-blk.centre[2])*(5.4/1000.0))
                if initT <= Tmin: initT=Tmin
            initSG=0.0 # initial gas saturation
            infvol=False # already given 1e50 volume
            pmx=grid.rocktype[rocktype].permeability[0]
        else:
            rocktype = 'main '
            initP=5e4+(997.0479*9.81*abs(hmax-blk.centre[2]))
            initSG=0.0
            initT=15.0+((np.abs(hmax-blk.centre[2])/100.0)*3.0)
            infvol=False
            if lay==geo.layerlist[-1]:
                rocktype='sourc'
            if (lpregion is not None and 'lp   ' in grid.rocktype.keys() and
                blk.centre[2] > lpregion[0][2] and 
                blk.centre[2] <= lpregion[1][2] and 
                blk.centre[0] > lpregion[0][0] and 
                blk.centre[0] <= lpregion[1][0]): # if in lp region
                rocktype='lp   '     
            if (hpregion is not None and 'hp   ' in grid.rocktype.keys() and
                blk.centre[2] > hpregion[0][2] and 
                blk.centre[2] <= hpregion[1][2] and 
                blk.centre[0] > hpregion[0][0] and 
                blk.centre[0] <= hpregion[1][0]): #if in hp region
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
            pmx=pmxcalc(blk,grid,hmax,rocktype,0.004,800.)      
        ptg.rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx,infvol=infvol)
    return grid
                
    
def pmxcalc(blk,grid,hmax,rock,Saar_lam=0.004,switch_depth=None):
    depth=hmax-blk.centre[2]
    if switch_depth is not None:
        if depth < switch_depth:    
            pmx=grid.rocktype[rock].permeability[0]*np.exp(-Saar_lam*depth) # calculating depth dependent permeability modifier
        else: pmx=grid.rocktype[rock].permeability[0]*np.exp(-Saar_lam*(switch_depth))*((depth/switch_depth)**-3.2)
    else:
        pmx=grid.rocktype[rock].permeability[0]*np.exp(-Saar_lam*depth)
    return pmx        

def heatgen(mod,geo,dat,grid,heat_flux):
    dat.clear_generators()
    f = open(mod+'/genertot.txt','w')
    f.write('Model = '+mod+'\n')
    allgens=[]
    cols=[col for col in geo.columnlist]
    f.write('Constant generation ='+str(heat_flux)+' J/s/m2\n')
    for col in cols:
        lay=geo.layerlist[-1] # bottom layer
        blkname=geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
        if grid.block[blkname].hotcell is not True:
            gxa=col.area*heat_flux
            gen=t2generator(name=' H'+col.name,block=blkname,type='HEAT',gx=gxa) # creat a generater oject with the heat generation rate of tflux - muliplication by column area important. 
            dat.add_generator(gen) # add generater to TOUGH2 input
            allgens.append(gxa)
    allgens=np.array(allgens)
    gensum=np.sum(allgens)
    f.write('Total generation in model = '+str(gensum)+' J/s\n')
    f.write('Total generation rate per m2 = '+str(gensum/geo.area)+' J/s/m2\n')
    f.close()

def simple_readres( modelname, savevtk=False, tough2_input=None, geom_data=None, results=None):
    """Easy reading of TOUGH2 results into pyTOUGH.  Option to write out vtk results files. Will return pyTOUGH results object if desired.
    
    Usage: 
    
    results=simple_readres(modelname, kwargs)
    
    modelname = directory name of model as string, e.g. "20150202_1"
    
    
    Optional kwargs:
    
    savevtk: flag to save vtk results files [default=False]
    
    tough2_input: provide the name of the tough2 input file [dafault="flow2.inp"]
    
    geom_data: provide name of pyTOUGH geometry data file ["grd.dat"]
    
    results: if results object already exists [None]
    
    """    
    mod=modelname # define modelname
    current_d=os.getcwd() # find current directory 
    print os.getcwd() 
    os.chdir(mod) # change to model directory
    print('Reading '+mod+' results from: ') #
    print os.getcwd()
    ## read output file
    t0=time.clock()
    if tough2_input is None:
       dat=t2data('flow2.inp') # tough2 input from input file
    else: dat=tough2_input

    if geom_data is None:
        geo=mulgrid('grd.dat') # geometry from geometry file
    else: geo=geom_data    

    if results is None:
        resfile='flow2.out'
#        resfile=gzip.open('flow2.out.gz','r').read()
        results=t2listing(resfile) # read output file
    else:
        results=results
    
    grid=dat.grid # define input grid
    t1=time.clock()
    t=t1-t0
    print 'time2read .out=',t
    
    ## create directory for results
    if not os.path.exists('results'):
        os.makedirs('results')
    os.chdir('results')
    
    ## save vtks...
    if savevtk:
       t0=time.clock()
       results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
       t1=time.clock()
       t=t1-t0
       print 'time2writevtks',t
    os.chdir(current_d)   
    return results