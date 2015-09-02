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
t0=tinit=time.clock()
plt.close('all')

read=True ########### I N P U T #########################
readgeo=True ########### I N P U T #########################
geo_fname='grd.dat'
readdat=True ########### I N P U T #########################
dat_fname='flow2.inp'
readresults=True ########### I N P U T #########################
results_fname='flow2.out'
readsat=True ########### I N P U T #########################
sat_fname='results/sat.pkl'


save=True ########### I N P U T #########################
savevtk=False ########### I N P U T #########################
batch_or_straight='st' ########### I N P U T #########################
modelorigin=(586034.886,1852660.465)


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

class ValidationError(Exception):
    def __init__(self, message, errors):
        # Call the base class constructor with the parameters it needs
        super(ValidationError, self).__init__(message)
        # Now for your custom code...
        self.errors = errors    

#%% setup
if batch_or_straight in anas+['b','ba','bat','batc']:
    batch=True
    param='wt' ########### I N P U T #########################
    main=True ########### I N P U T #########################
else:
    batch=False
    mod='20150806_2_var1' ########### I N P U T #########################


#%%###########################################################################
############################# STRAIGHT MODE ###############################
if not batch: 
    t0=time.clock()
    print 'running in straight mode (',batch_or_straight,')'
    print 'model=',mod
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/'+mod.split('_var')[0]+'/'+mod)
    if read:
        if readgeo is True: 
            print 'Reading geometry from '+ geo_fname
            geo=mulgrid(geo_fname)
            geo.radial=False # set flag that geometry is not yet radial.
        if readdat is True: 
            print 'Reading input data from '+ dat_fname
            dat=t2data(dat_fname) 
        if readsat is True: 
            if os.path.isfile(sat_fname):
                print 'Reading saturation data from '+ sat_fname
                sat=ptg.load_obj(sat_fname)
            else: 
                print('CANT READ SAT FILE......Continuing with sat={}')
                sat={}                
        if readresults is True:
            if sat=={}:
                print 'Reading results from '+ results_fname
                results=t2listing(results_fname)
            else:
                results=[]
            

                
                
    t1=time.clock()        
    print 'time to read=',(t1-t0)
    width=geo.bounds[1][1]-geo.bounds[0][1] #10.0
    ## define well locations
#    stations=np.genfromtxt('../dev_files/Steffie_station_locs.txt', delimiter=',', dtype=None, skiprows=1, usecols=(0,1) ,names='x,y')
#    station_dists=[np.sqrt(((modelorigin[0]-x)**2 + (modelorigin[1]-y)**2)) for x,y in zip(stations['x'],stations['y'])]
#    station_dists.sort()
#    wellx=station_dists
    wellx=wellx=[50.0,500.0,1100.0,1750,2650,3200]#[50.0,500.0,1750,2200,2650]##5.0,50.0,500,1750,2200,2650#######################################################
    #wellx=[50.0]################################################################
    welly=[width/2.0]*len(wellx) # make y coords same length as x coords
    wells=np.hstack((np.transpose([wellx]),np.transpose([welly])))
    t2=time.clock()
    t=t2-t0
    print 'time2setup=',t
    results,sat,wellsatblk=ptg.readres(mod,wells,save=save,savevtk=savevtk,results=results,sat=sat,tough2_input=dat, geom_data=geo, maxtime=100)
      
       
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
        if read:
            if readgeo is True: 
                print 'Reading geometry from '+ geo_fname
                geo=mulgrid(geo_fname)
            if readdat is True: 
                print 'Reading input data from '+ dat_fname
                dat=t2data(dat_fname) 
            if readresults is True: 
                print 'Reading results from '+ results_fname
                results=t2listing(results_fname)
        ptg.readres(mod,wells,save=save,savevtk=savevtk,results=results,tough2_input=dat, geom_data=geo,fall=fall)
    fall.close()
print 'time to run =', time.clock()-tinit
