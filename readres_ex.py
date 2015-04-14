# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 18:41:58 2014
Script for reading results calculating and plotting gravity

@author: glbjch
"""

import os
os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough')
import pytoughgrav as ptg
import numpy as np
import time
import matplotlib.pyplot as plt
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *

#%% set up
t0=time.clock()
plt.close('all')

read=False ########### I N P U T #########################
save=True ########### I N P U T #########################
savevtk=True ########### I N P U T #########################
batch_or_straight='st' ########### I N P U T #########################
modelorigin=(586034.886,1852660.465)
width=10.0

#%% functions
def anagrams(word):
    """ Generate all of the anagrams of a word. """ 
    if len(word) < 2:
        yield word
    else:
        for i, letter in enumerate(word):
            if not letter in word[:i]: #avoid duplicating earlier words
                for j in anagrams(word[:i]+word[i+1:]):
                    yield j+letter 

anas=[]
if __name__ == "__main__":
    for i in anagrams("batch"):
        anas=anas+[i]
        

def readfiles(read,geo,dat,results):
    t0=time.clock()
    if read:
        geo=mulgrid('grd.dat') ########### I N P U T #########################
        dat=t2data('flow2.inp') ########### I N P U T #########################
        results=t2listing('flow2.out')
    else:
        try: geo 
        except NameError: geo=None
        try: dat
        except NameError: dat=None
        try: results
        except NameError: results=None
        #grid=dat.grid # define input grid
    return geo,dat,results
    print 'time to read=',(time.clock()-t0)

#%% setup
if batch_or_straight in anas+['b','ba','bat','batc']:
    batch=True
    param='wt' ########### I N P U T #########################
    main=True ########### I N P U T #########################
else:
    batch=False
    mod='20150327_1_var' ########### I N P U T #########################


    

#%%###########################################################################
############################# STRAIGHT MODE ###############################
if not batch: 
    t0=time.clock()
    print 'running in straight mode (',batch_or_straight,')'
    print 'model=',mod
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Steffi_GRAV/'+mod)    
    geo,dat,results=readfiles(read,None,None,results)#,geo,dat,results)
    ## define well locations
    stations=np.genfromtxt('../dev_files/Steffie_station_locs.txt', delimiter=',', dtype=None, skiprows=1, usecols=(0,1) ,names='x,y')
    station_dists=[np.sqrt(((modelorigin[0]-x)**2 + (modelorigin[1]-y)**2)) for x,y in zip(stations['x'],stations['y'])]
    station_dists.sort()
    wellx=station_dists
    #wellx=[50.0,500.0,1750,2200,2650]##5.0,50.0,500,1750,2200,2650#######################################################
    #wellx=[50.0]################################################################
    welly=[width/2.0]*len(wellx) # make y coords same length as x coords
    wells=np.hstack((np.transpose([wellx]),np.transpose([welly])))
    t1=time.clock()
    t=t1-t0
    print 'time2setup=',t
    results=ptg.readres(mod,wells,save=save,savevtk=savevtk,results=results,tough2_input=dat, geom_data=geo)

#%%################################################################
######################### BATCH MODE ############################
if batch:
    print 'running in batch mode (',batch_or_straight,')'
    print 'parameter=',param
    if main:
        print 'working on main model results'
    else:
        print 'working on warm up model results'
    #os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
    if save:
       fall=open('ampslist.txt','w')
       fall.write('modelname \t massmin \t massmax \t massamp\t gravmin \t gravmax \t gravamp \n')
    wellx=[0.5]
    welly=[0.5]*len(wellx)
    wells=hstack((np.transpose([wellx]),np.transpose([welly])))
    for mod in [mod for mod in os.listdir('../'+param) if os.path.isdir(os.path.join('../'+param,mod))]:
        #os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
        os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
        if main:
           os.chdir(mod+'/main')
        else:
           os.chdir(mod)
        geo,dat,results=readfiles(read)
        ptg.readres(mod,wells,save=save,savevtk=savevtk,results=results,tough2_input=dat, geom_data=geo,fall=fall)
    fall.close()
print 'time to run =', time.clock()-t0
       