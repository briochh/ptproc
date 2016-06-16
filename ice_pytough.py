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
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import copy
import matplotlib.colors

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

def icegrid(geo,dat,rocks,boundcol,eos=3,lpregion=None,hpregion=None,heatsource=None,satelev=0.0,atmosP=1.013e5,pmx_lamda=0.004, glacier_limit=2500., infax=False, topsurf=None,radial=False):
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
        if topsurf is not None:
            ind=topsurf[:,0].tolist().index(col.centre[0])
            hmax=topsurf[ind,1]
        else: hmax=col.surface
        if blk in atmos:
            rocktype='top  ' # assign rocktype "top  "
            initP=atmosP*spy.power(1.-(col.surface*2.25577e-5),5.25588)  # initial presure condition - MAY NOT BE APPROPRIATE - WHAT IS THE PRESSURE UNDER THICK GLACIER AT 5000 m amsl??         
            Tmin=2.
            #if blk.centre[0] <= 250.:
                #initT=50.#Tmin + (hmax-blk.centre[2])*(12.5/100.)
            if 350.< blk.centre[0] <= glacier_limit: 
                initT=Tmin#+(hmax-blk.centre[2])*(12.5/100.) # initial temperature - TOUGH2 doesn't seem to like < 1.0 C
            else:
                initT = 25.8 - (hmax*(5.4/1000.)) #+ (hmax-blk.centre[2])*(12.5/100.) # 15.+((2000.-blk.centre[2])*(5.4/1000.0))
            if initT <= Tmin: initT=Tmin
            #initT=25.8 - (hmax*(5.4/1000.)) + (hmax-blk.centre[2]*15./100.) 150.+(-0.0225*blk.centre[0])
            if (hpregion is not None and 'hp   ' in grid.rocktype.keys()):
                for rt,hpr in hpregion.iteritems():
                    if (blk.centre[2] > hpr[0][2] and 
                    blk.centre[2] <= hpr[1][2] and 
                    blk.centre[0] > hpr[0][0] and 
                    blk.centre[0] <= hpr[1][0]): #if in hp region
                        rocktype=rt#'hp   ' # this allows a different pmx for atmos above highperm
            initSG=0.999 # initial gas saturation
            infvol=False # already given 1e50 volume
#            if topsurf is not None:
#                ind=topsurf[:,0].tolist().index(col.centre[0])
#                hmax=topsurf[ind,1]
#                pmx=pmxcalc(blk,grid,hmax,rocktype,0.004,800.)
#            else:
            pmx=grid.rocktype[rocktype].permeability[0]
            rocktype='top  ' # resets to rocktype "top  "
        else:
            rocktype = 'main '
            initP=(atmosP*spy.power(1.-(col.surface*2.25577e-5),5.25588))+(997.0479*9.81*abs(col.surface-blk.centre[2]))
            if blk.centre[2]<4800:
                initSG=0.0
            else:
                initSG=10.3
            if blk.centre[0]<10000:
                initT=Tmin + ((np.abs(hmax-blk.centre[2])/100.0)*12.5)
            else:            
                initT=Tmin + 10.0 + ((np.abs(hmax-blk.centre[2])/100.0)*3.0)
            if initT > 350.: initT=350.
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
                for rt,hpr in hpregion.iteritems():
                    if (blk.centre[2] > hpr[0][2] and 
                    blk.centre[2] <= hpr[1][2] and 
                    blk.centre[0] > hpr[0][0] and 
                    blk.centre[0] <= hpr[1][0]): #if in hp region
                        rocktype=rt#'hp   '
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
            if not radial and col in boundcol:
                print "inf vol boundary cell " + blk.name
                infvol=True
                initSG=0.0
                rocktype='bound'
            pmx=pmxcalc(blk,grid,hmax,rocktype,0.004,800.)      
        ptg.rockandincon(blk,grid,dat,rocktype,initP,initSG,initT,pmx,eos=eos,infvol=infvol)
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

def heatgen(mod,geo,dat,grid,heat_flux,function=None, inject=None, inject2=None):
    f = open(mod+'/Heat_genertot.txt','w')
    f.write('Model = '+mod+'\n')
    allgens=[]
    allinject=[]
    cols=[col for col in geo.columnlist]
    if function is not None:
        func=function['type']
        if func=='exp':
            p1=function['points'][0]
            p2=function['points'][1]
            b=spy.log(p1[1]/p2[1])/(p1[0]-p2[0])
            a=p1[1]/spy.e**(p1[0]*b)
            f.write('exponential spatial generation Q=ae^bx '+str(a)+'e^('+str(b)+'x) J/s/m2\n')
        elif func=='log':
            p1=np.array(function['points'][0])
            p2=np.array(function['points'][1])
            a=(p1[1]-p2[1])/(spy.log(p1[0]/p2[0]))            
            b=spy.e**((p2[1]*spy.log(p1[0])-p1[1]*spy.log(p2[0]))/(p1[1]-p2[1]))
            f.write('logarithmic spatial generation Q=a*ln(bx) '+str(a)+'ln('+str(b)+'x) J/s/m2\n')
    else: 
        func='Constant'
        f.write('Constant generation ='+str(heat_flux)+' J/s/m2\n')
    for col in cols:
        lay=geo.layerlist[-1] # bottom layer
        blkname=geo.block_name(lay.name,col.name) # get block name for the bottom layer of this column
        if grid.block[blkname].hotcell is not True:
            if func=='exp':
                heat_flux=a*spy.e**(b*grid.block[blkname].centre[0])
            elif func=='log':
                #if grid.block[blkname].centre[0] < 250:
                #    heat_flux=p1[1]
                #else:
                    heat_flux=a*spy.log(b*(grid.block[blkname].centre[0]))
            #gxa=[0.0]+[col.area*heat_flux]*2
            gxa=col.area*heat_flux
            #times=[0.]+[1000*365.25*3600*24]+[1.0e6*365.25*3600*24]
            #numt=len(times)
            gen=t2generator(name=' H'+col.name,block=blkname,type='HEAT',gx=gxa, ex=None,hg=None,fg=None)#, rate=gxa, time=times, ltab=numt) # creat a generater oject with the heat generation rate of tflux - muliplication by column area important. 
            dat.add_generator(gen) # add generater to TOUGH2 input
            allgens.append(gxa)
            if inject is not None:
                if grid.block[blkname].centre[0] < inject[0]:
                    ixa=col.area*inject[1]
                    gen=t2generator(name=' i'+col.name,block=blkname,type='COM1',gx=ixa, ex=inject[2]) # create a generater oject with the heat generation rate of tflux - muliplication by column area important. 
                    dat.add_generator(gen)
                    allinject.append(ixa)
    if inject2 is not None:
        lays=[lay for lay in geo.layerlist if inject2[0] < lay.centre < inject2[1]]
        for lay in lays:
            col=geo.columnlist[0] # first column
            blkname=geo.block_name(lay.name,col.name)
            ixa=lay.get_thickness()*2*np.pi*col.get_side_lengths()[1]*inject2[2]
            gen=t2generator(name='i'+lay.name+'01',block=blkname,type='COM1',gx=ixa, ex=inject2[3]) # create a generater oject with the heat generation rate of tflux - muliplication by column area important. 
            dat.add_generator(gen)
            allinject.append(ixa)
    allgens=np.array(allgens)
    gensum=np.sum(allgens)    
    allinject=np.array(allinject)
    injectsum=np.sum(allinject)
    f.write('Total generation in model = '+str(gensum)+' J/s\n')
    f.write('Total generation rate per m2 = '+str(gensum/geo.area)+' J/s/m2\n')
    if inject is not None:
        f.write('Total injection in model = '+str(injectsum)+' kg/s\n')
        f.write('injection rate per m2 = '+str(inject[1])+' kg/s/m2\n')
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
    return geo,grid,dat,results
    
def icepost( modelname, save=False, radial=True, savevtk=False, geom_data=None, tough2_input=None, results=None, times={}, fall=None,flows={'FLOH':{},'FLOF':{}}, logt=False):
    """ Function to calculated surface heat flow from pytough results
    """
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
        print('No Results files (flow.out) passed. please read flow2.out')
    #results=results   
    width=geo.bounds[1][1]-geo.bounds[0][1]   
    yrsec=365.25*3600*24
    mod=modelname
    glaclim=[250.,2500.]
    if radial and not geo.radial:    
        ptg.makeradial(geo,None,width)
    # find atmosphere blocks
    grid=dat.grid # define input grid    
    grid2=t2grid().fromgeo(geo) # grid saved in flow2.inp does not contain atmos information required.
    atmos=[]
    atmosconn=[] 
    os.chdir('results')
    for flow in flows:
        t0=time.clock()
        if flows[flow]=={}:
            if atmos==[]:
                atmos=grid2.atmosphere_blocks
                # find atmosphere connections
                atmosconn=[cell.connection_name for cell in atmos] # returns a list of sets - !!!! sets are not ordered so cannot be indexed !!! 
            flows[flow],times = surfaceflow(atmos,atmosconn,results=results,grid=grid2,flow=flow)
        tq=times
        tscale=tq/yrsec
        if logt:
            tscale=np.log10(tscale) 

        # produce flow matrix (qts IN KG/S/M2 OR W/MS)
        X,Area,qts = prod_flowmatrix(flows[flow].values())       
                
        delref=np.array(qts.T[0]) # reference outflow at t0 kg/s/m2 ~ mm/s # mask where flow is negative (i.e. in to the model)
        craterflowts,_ =convert2rate(qts.T,tscale,[0,250],Area,X)
        glacflowts,_ =convert2rate(qts.T,tscale,glaclim,Area,X)
        
        
        if flow=='FLOH' or flow=='FHEAT': # calculate meltrate etc.
            unit='W'
            # compute meltrate in surface cells from heat out flow
            # some meltrate plots performed in this one
            calcmeltrate(mod,qts,delref,X,tscale,Area,
                         glaclim=[250.,2500.],save=save) 
        else:
            unit='kg/s'
            
        plotflows(mod,flow,qts,X,tscale,Area,unit,delref,glaclim,save=save)
#        # plottings

        if save:
            if os.path.isfile(flow+'.pkl'):
                print(flow+' flow alreay pickled')
            else:
                ptg.save_obj(flows[flow],flow+'.pkl')  
            if os.path.isfile('time.pkl'):
                print('time alreay pickled')
            else:
                ptg.save_obj(tq,'time.pkl')
        t1=time.clock()   
        print 'time for flow=',(t1-t0)  
            
    if savevtk and os.path.isfile(mod+'_output.vtk') is False:
        results.write_vtk(geo,mod+'_output.vtk',grid=grid,flows=True)
        #os.chdir('..')
    
def plotflows(mod,flow,qts,X,tscale,Area,unit,ssqts,glaclim=[250.,2500.],save=True):
    if ssqts is not None:
        dodelta=True
    else: 
        dodelta=False
    qin=qts.clip(max=0.)# just where flow is where flow is negative (i.e.in to the model)
    qout=qts.clip(min=0.)# just where flow is where flow is positive (i.e.out of the model)
    if dodelta:
        ssqout=ssqts.clip(min=0)
        ssqin=ssqts.clip(max=0)
        delqout=np.subtract(qout.T,ssqout) # kg/s/m2 ~ mm/s
        delqin=-np.subtract(qin.T,ssqin) 
    craterflow=qts[X<glaclim[0]] # W/m2 or kg/s/m2
    glacflow=qts[(X<glaclim[1]) & (X>glaclim[0])]   
    craterflowts,_ =convert2rate(qts.T,tscale,[0,250],Area,X)
    glacflowts,glac_contactArea =convert2rate(qts.T,tscale,glaclim,Area,X)

    # plottings
    ## a quick plot of flows into atmosphere at X and time.
    plt.figure() # FLOW
    plt.pcolormesh(X,tscale,qts.T, rasterized=True,cmap='rainbow',norm=matplotlib.colors.LogNorm(vmin=1e-2,vmax=100)) #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.1e')
    #cbar.set_label(flow + r' out of the model ('+ unit + r'/m$^{2}$)')
    cbar.set_label('Surface heat flux ('+ unit + r'm$^{-2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.ylim(tscale.min(),5000)
    #plt.title('Flow ('+flow+') out of the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    plt.tight_layout()
    if save:
        plt.savefig(mod+'_'+flow+'_.pdf',dpi=400)
    else:
        plt.savefig(mod+'_'+flow+'_log_.pdf',dpi=400)
    
    #qout=np.ma.masked_array(qts,[qts<0]) # mask where flow is negative (i.e. in to the model)
    #qin=np.ma.masked_array(qts,[qts>0]) # mask where flow is positive (i.e. out of the model)

    if dodelta:
        plt.figure()
        plt.pcolormesh(X,tscale,np.ma.masked_array(delqout,[qts.T<0]), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
        cbar=plt.colorbar(format='%.1e')
        cbar.set_label('Change in '+ flow + r' out of the model ('+ unit + r'/m$^{2}$)')
        #plt.xlim(0,2500)
        plt.ylim(tscale.min(),tscale.max())
        plt.title('Change in flow ('+flow+') out of the model')
        plt.xlabel('Distance from centre axis (m)')
        plt.ylabel('Time (yrs)')
        plt.tight_layout()
        if save:
            plt.savefig(mod+'_'+flow+'_delout_.pdf',dpi=400)
        
        plt.figure()
        plt.pcolormesh(X,tscale,np.ma.masked_array(delqin,[qts.T>0]), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
        cbar=plt.colorbar(format='%.1e')
        cbar.set_label('Change in '+ flow + r' in to the model ('+ unit + r'/m$^{2}$)')
        #plt.xlim(0,2500)
        plt.ylim(tscale.min(),tscale.max())
        plt.title('Change in flow ('+flow+') in to the model')
        plt.xlabel('Distance from centre axis (m)')
        plt.ylabel('Time (yrs)')
        plt.tight_layout()
        if save:
            plt.savefig(mod+'_'+flow+'_delin_.pdf',dpi=400)        
    #        plt.figure()


    crmesh=plt.figure()
    plt.pcolormesh(X[X<glaclim[0]],tscale,np.ma.masked_array(craterflow,[craterflow<0]).T, rasterized=True,cmap='rainbow') # mm/s
    cbar=plt.colorbar(format='%.1e')
    cbar.set_label('Flow ('+flow+') at crater')
    plt.xlim((0,glaclim[0]))
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Crater flow ('+flow+')')
    
    crts=plt.figure()
    plt.plot(tscale,craterflowts)
    #plt.xlim(0.08, 100)
    #plt.ylim(craterflowts.min(),craterflowts.max()) #!!!!!!!!
    plt.xlabel('Time (yrs)')
    plt.ylabel('Crater Flow ('+unit+')')
    plt.title('Crater flow ('+flow+')')
    if save:
        crts.savefig(mod+'_'+flow+'_crater_ts.pdf',dpi=400)
        crmesh.savefig(mod+'_'+flow+'_craterflow.pdf',dpi=400) 
        np.savetxt(mod+'_'+flow+'_CraterFlow.txt',np.vstack(([tscale],[craterflowts])).T)
        
    vmax=5.0 # allow quick change colour limit
    #vmax=320
    fmat='%.1'  # allow quick change of number format
    
    glacmesh=plt.figure()
    plt.pcolormesh(X[(X<glaclim[1]) & (X>glaclim[0])],tscale,np.ma.masked_array(glacflow,[glacflow<0]).T, rasterized=True,cmap='rainbow',vmin=0.0,vmax=vmax) # mm/s
    cbar=plt.colorbar(format=fmat+'f')
    cbar.set_label('Flow ('+flow+') into galcier ('+ unit + r'm$^{-2}$)')
    plt.xlim((0,2500))
    plt.ylim(0,1000) # allow focus to particular timescale
    plt.ylabel('Time (yrs)')
    plt.xlabel('Distance x (m)')
    #plt.ylim(tscale.min(),tscale.max())
    plt.title('Glacier flow ('+flow+')')
    
    
    #glacts=plt.figure()
    glacts,ax1=plt.subplots()
    plot_ax1, = ax1.plot(tscale,glacflowts)
    ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
    #plt.ylim(craterflowts.min(),craterflowts.max()) #!!!!!!!!
    ax1.set_xlabel('Time (yrs)')
    ax1.set_ylabel('Glacier Flow ('+unit+')')
    ax1.tick_params(axis='y', colors=plot_ax1.get_color())
    ax1.yaxis.label.set_color(plot_ax1.get_color())
    ax1.spines['left'].set_color(plot_ax1.get_color())
    ax1.yaxis.get_offset_text().set_color(plot_ax1.get_color())
    ax2=plt.twinx(ax1)
    plot_ax2, = ax2.plot(tscale,glacflowts/glac_contactArea,'g')
    ax1.set_xlim(0, 1000)
    ax1.set_ylim(0, 7.0e5)
    ax2.set_ylim(0, 3.5e-2)
    #ax2.set_ylim(0, 250)
    ax2.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,3))
    ax2.tick_params(axis='y', colors=plot_ax2.get_color())
    ax2.yaxis.label.set_color(plot_ax2.get_color())
    ax2.spines['right'].set_color(plot_ax2.get_color())  
    ax2.yaxis.get_offset_text().set_color(plot_ax2.get_color())
    ax2.set_ylabel(r'Flux to glacier ('+ unit + r'm$^{-2}$)')
    plt.title('Glacier flow ('+flow+')')
    if save:
        glacts.savefig(mod+'_'+flow+'_glacier_ts.pdf',dpi=400,bbox_inches='tight')
        glacmesh.savefig(mod+'_'+flow+'_glacierflow.pdf',dpi=400,bbox_inches='tight') 
        np.savetxt(mod+'_'+flow+'_GlacFlow.txt',np.vstack(([tscale],[glacflowts])).T)
        
    glacmesh200,ax200=plt.subplots()
    im=ax200.pcolormesh(X[(X<glaclim[1]) & (X>glaclim[0])],
                     tscale,np.ma.masked_array(glacflow,[glacflow<0]).T,
                     rasterized=True,cmap='rainbow',
                     norm=matplotlib.colors.LogNorm(vmin=1e-2,vmax=vmax)) # mm/s

    cbar=glacmesh200.colorbar(im, ticks=LogLocator(subs=range(10)))#,format=fmat+'e')
    cbar.set_label('Flow ('+flow+') into galcier ('+ unit + r'm$^{-2}$)')
    ax200.set_xlim((0,2500))
    ax200.set_ylim(0,200)
    ax200.set_ylabel('Time (yrs)')
    ax200.set_xlabel('Distance x (m)')
    #plt.ylim(tscale.min(),tscale.max())
    ax200.set_title('Glacier flow ('+flow+')')
    
    glacts200,ax1=plt.subplots()
    plot_ax1, = ax1.plot(tscale,glacflowts)
    ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
    #plt.ylim(craterflowts.min(),craterflowts.max()) #!!!!!!!!
    ax1.set_xlabel('Time (yrs)')
    ax1.set_ylabel('Glacier Flow ('+unit+')')
    ax1.tick_params(axis='y', colors=plot_ax1.get_color())
    ax1.yaxis.label.set_color(plot_ax1.get_color())
    ax1.spines['left'].set_color(plot_ax1.get_color())
    ax1.yaxis.get_offset_text().set_color(plot_ax1.get_color())
    ax2=plt.twinx(ax1)
    plot_ax2, = ax2.plot(tscale,glacflowts/glac_contactArea,'g')
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 5.0e5)
    ax2.set_ylim(0, 2.5e-2)
    #ax2.set_ylim(0, 250)
    ax2.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,3))
    ax2.tick_params(axis='y', colors=plot_ax2.get_color())
    ax2.yaxis.label.set_color(plot_ax2.get_color())
    ax2.spines['right'].set_color(plot_ax2.get_color())  
    ax2.yaxis.get_offset_text().set_color(plot_ax2.get_color())
    ax2.set_ylabel(r'Flux to glacier ('+ unit + r'm$^{-2}$)')
    plt.title('Glacier flow ('+flow+')')
    if save:
        glacts200.savefig(mod+'_'+flow+'_glacier_ts200.pdf',dpi=400,bbox_inches='tight')
        glacmesh200.savefig(mod+'_'+flow+'_glacierflow200.pdf',dpi=400,bbox_inches='tight') 


def prod_flowmatrix(values,prelen=0):
        X=[]
        qts=[]
        Area=[]
        for x,a,q in values:
            X.append(x)
            qts.append(q[prelen:])
            Area.append(a)
            #ssqts.append(q[0])
        inds=np.array(X).argsort()
        X=np.array(X)[inds]
        Area=np.array(Area)[inds]
        qts=np.array(qts)[inds] # J/s/m2 or kg/s/m2 
        #ssqts=np.array(ssqts)[inds]
        #qin=np.ma.masked_array(qts,[qts>0]) # mask where flow is positive (i.e. out of the model)
        #qout=np.ma.masked_array(qts,[qts<0]) # mask where flow is negative (i.e. in to the model)
        return X,Area,qts#,ssqts#,qin,qout
    
def convert2rate(matrix,times,Xlimits,Areas,Xlocs):
    timeseries=np.zeros(len(times))
    TotArea=0
    i=0
    for t in times:
            for x,A,r in zip(Xlocs,Areas,matrix[i]):
                if r > 0 and x < Xlimits[1] and x > Xlimits[0]:
                    timeseries[i]= timeseries[i] + (r*A) # kg/s or W
                    if i==0: TotArea=TotArea+A
            i+=1
    return timeseries,TotArea

def calcmeltrate(mod,qts,ssqts,X,tscale,Area,glaclim=[250,2500],save=True,logt=False,plot=True):
        #unit='W'
        if ssqts is not None:
            dodelta=True
        else: 
            dodelta=False
        qout=qts.clip(min=0)
        yrsec=365.25*24.*3600.
        meltratematrix=(qout.T/3.34E5) # kg/s/m2
        #meltratematrix[meltratematrix<0]=0 # kg/s/m2
        if dodelta: 
            ssmelt=(ssqts.clip(min=0.)/3.34E5)
            # change in meltrate
            deltameltrate=np.subtract(meltratematrix,ssmelt) # kg/s/m2 ~ mm/s
            deltaglacmeltrate=deltameltrate.T[(X>glaclim[0]) & (X<glaclim[1])].T # kg/s/m2 ~ mm/s
        #meltrate just within glacier            
        glacmeltrate=meltratematrix.T[(X>glaclim[0]) & (X<glaclim[1])].T  # kg/s/m2 ~ mm/s
        #change in meltrate within glacier 
        meltrate,glacArea=convert2rate(meltratematrix,tscale,glaclim,Area,X)
        
#        meltrate=np.zeros(len(tscale))
#        glacArea=0
#        i=0
#        for t in tscale:
#            for x,A,r in zip(X,Area,meltratematrix[i]):
#                if r > 0 and x < glaclim[1] and x > glaclim[0]:
#                    meltrate[i]= meltrate[i] + (r*A) # kg/s
#                    glacArea=glacArea+A
#            i=i+1  
        meltrate_mmpyr= (meltrate*yrsec)/glacArea#((np.pi*(glaclim[1]**2))-(np.pi*(glaclim[0]**2))) # kg/yr/m2 ~ mm/yr 
        meltvol=[rate*((tscale[i+1]-tscale[i])*yrsec)/1000. for rate,i in zip(meltrate,range(len(tscale)-2))]
        meltvolcum=np.cumsum(meltvol)
        basemeltlvol=meltrate[0]*(tscale*yrsec)/1000.+meltvolcum[0]
        difmeltvolcum=np.subtract(meltvolcum,basemeltlvol[0:-2])
        if difmeltvolcum.max() >=1.0e5:
            ind=np.where(difmeltvolcum>=1e5)[0][0]
            print "time when melt volume is 1e5 m3 =",tscale[ind]
        else:
            print "Melt volume never exceeds 1e5 m3"
        ####  plotsss
        if plot:
            plt.figure()
            plt.pcolormesh(X,tscale,meltratematrix, rasterized=True,cmap='rainbow') # mm/s
            cbar=plt.colorbar(format='%.1e')
            cbar.set_label(r'Melt rate (kg/s/m$^{2}$)')
            #plt.xlim((0,2500))
            plt.ylim(tscale.min(),tscale.max())
            plt.xlabel('Distance from centre axis (m)')
            plt.ylabel('Time (yrs)')    
            plt.title('Melt rate')
            plt.tight_layout()
            if save:
                plt.savefig(mod+'_meltrate_.pdf',dpi=400)  
        
            plt.figure()
            plt.pcolormesh(X[(X>glaclim[0]) & (X<glaclim[1])],tscale,glacmeltrate, rasterized=True,cmap='rainbow')#, vmax=3.2e-4) # mm/s
            cbar=plt.colorbar(format='%.1e')
            cbar.set_label(r'Glacial melt rate (kg/s/m$^{2}$)')
            #plt.xlim((0,250))
            plt.ylim(tscale.min(),tscale.max())
            plt.xlabel('Distance from centre axis (m)')
            plt.ylabel('Time (yrs)')    
            plt.title('Melt rate')
            plt.tight_layout()
            if save:
                plt.savefig(mod+'_glac_meltrate_.pdf',dpi=400)  
           
            if dodelta:
                plt.figure()
                plt.pcolormesh(X,tscale,deltameltrate, rasterized=True,cmap='rainbow') # mm/s
                cbar=plt.colorbar(format='%.1e')
                cbar.set_label(r'Change in melt rate (kg/s/m$^{2}$)')
                #plt.xlim((0,2500))
                plt.ylim(tscale.min(),tscale.max())
                plt.xlabel('Distance from centre axis (m)')
                plt.ylabel('Time (yrs)') 
                plt.title('Change in melt rate')
                plt.tight_layout()
                if save:
                    plt.savefig(mod+'_meltrate_delta_.pdf',dpi=400)  
            
                plt.figure()
                plt.pcolormesh(X[(X>glaclim[0]) & (X<glaclim[1])],tscale,deltaglacmeltrate, rasterized=True,cmap='rainbow') # mm/s
                cbar=plt.colorbar(format='%.1e')
                cbar.set_label(r'Change in melt rate (kg/s/m$^{2}$)')
                plt.xlim((0,2500))
                plt.ylim(tscale.min(),tscale.max())
                plt.title('Change in glacial melt rate')
                plt.tight_layout()
                if save:
                    plt.savefig(mod+'_glacmeltrate_delta_.pdf',dpi=400) 
    
        #############
        

           
    #        plt.figure()
    #        plt.plot(tscale,meltrate_mmpyr)
    #        #plt.xlim(0, 30000)
    #        plt.xlabel('Time (yrs)')
    #        plt.ylabel('Average melt rate at glacial base (mm/yr)')
    #        plt.title('Average basal meltrate')
    #        plt.tight_layout()
    #        if save:
    #            plt.savefig(mod+'_basalmelt.pdf')
                
            plt.figure()
            plt.semilogy(tscale[0:-2],meltvolcum)
            plt.plot(tscale,basemeltlvol)
            plt.plot(tscale[0:-2],difmeltvolcum)
            plt.xlim(0, 1000)
            plt.ylim(0,1e7)
            plt.xlabel('Time (yrs)')
            plt.ylabel('cumulative melt volume')
            #plt.title('Average basal meltrate')
            plt.tight_layout()
            if save:
                plt.savefig(mod+'_basalmelt.pdf')
                
            plt.figure()
            plt.plot(tscale,3600.*meltrate/20.)
            plt.xlim(0, 25)
            #plt.ylim(0,4.5e7) #!!!!!!!!
            plt.xlabel('Time (yrs)')
            plt.ylabel('Rate of mass loss from glacier (m3/hour)')
            plt.title('rate of mass loss')
            if save:
                plt.savefig(mod+'_melt_ts.pdf')
    
        if save:
            np.savetxt("melt_time.txt", np.vstack(([tscale], [meltrate],[meltrate_mmpyr])).T)
        #return

        
def surfaceflow(atmoscells,atmosconns,results=None,grid=None,flow='FLOH'):
#    totq=glac_toq=np.zeros(results.num_times) # arrays for adding total flows at each timstep 
#    X=np.array([]).reshape(0,1) # array for storing the x dimension of each connection 
#    Area=np.array([]).reshape(0,1)
#    qts=[] # empty list for storing all flows
    flowDict={} 
    sels=[('c',conn,flow) for atmosset in atmosconns for conn in atmosset]
    Q=results.history(sels)
    tq=Q[0][0]
    for sel,q in zip(sels,Q):
        conn=sel[1]
        cell1,cell2=conn
        # tq is time, q is flow
        #(tq,q) = results.history([('c', conn, flow)]) # here just pull out flow of HEAT (could do more)
        # flow from cell1 to cell2 is negative we want to make it uniform so that flow into the atmosphere is always negative
        if grid.block[cell2] in atmoscells:
            q=-q[1]/(grid.connection[conn].area) # divide by conneciton area to get flux in /m2
            cen=grid.block[cell1].centre # 
        elif grid.block[cell1] in atmoscells:
            q=q[1]/(grid.connection[conn].area) # if atmos cell is cell1 the sign of the flow needs to be switched
            cen=grid.block[cell2].centre
        flowDict[conn]=cen[0],grid.connection[conn].area,q # dictionary of flow timseries for each connection

        # compile all the flow time series (for each connection)        
#            if qts == []:
#                qts=[q] # the first time through need to initialise our list
#            else:
#                qts=np.concatenate((qts,[q]),axis=0) # add the current flow time series to our compilations
#            X=np.append(X,cen[0]) # also create array of connection distance
#            Area= np.append(Area,grid2.connection[conn].area)
#            totq=np.add(totq,q) # timeseries of total heat flux into atmos 
#            # can look at total heat flux that is within a certain X limit (e.g. within the range of a glacier)        
#            if cen[0] <= 2500.0:
#                glac_toq=np.add(glac_toq,q) 
        #atxs[conn]=cen[0]
    return flowDict,tq
    
def makewt(model,geo=None,grid=None,dat=None,results=None):
    t0=tinit=time.clock()
    mod=model # define model name
    read=True ########### I N P U T #########################
    readgeo=True ########### I N P U T #########################
    geo_fname='grd.dat'
    readdat=True ########### I N P U T #########################
    dat_fname='flow2.inp'
    readresults=True ########### I N P U T #########################
    results_fname='flow2.out'
    
    print 'model=',mod
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+mod)  
    
    if read:
        if geo is None: 
            print 'Reading geometry from '+ geo_fname
            geo=mulgrid(geo_fname)
            geo.radial=False # set flag that geometry is not yet radial.
        if dat is None: 
            print 'Reading input data from '+ dat_fname
            dat=t2data(dat_fname) 
        if results is None:
            print 'Reading results from '+ results_fname
            results=t2listing(results_fname)
    t1=time.clock()        
    print 'time to read=',(t1-t0)  
    
    if grid is None:
        grid=dat.grid # define input grid    
    xz=[]
    results.last()
    sat=results.element['SL']
    for col in geo.columnlist:
        colblklist= [block for block in grid.blocklist if block.centre[0]==col.centre[0] and sat[geo.block_name_index[block.name]]==1.0]
        wtelev=np.max([geo.block_surface(geo.layer[geo.layer_name(str(blk))],col) for blk in colblklist])
        xz.append([col.centre[0], wtelev])
    xz.append([xz[-1][0]+xz[0][0],xz[-1][1]])
    xz.insert(0,[xz[0][0]-xz[0][0],xz[0][1]])

    xz=np.array(xz)

    np.savetxt('surface.txt',xz,header='/Distance\tHeight',fmt='%g',delimiter='\t')
    os.chdir('..')