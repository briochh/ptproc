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
import time
import ice_pytough as ipt

#%% set up
t0=tinit=time.clock()
plt.close('all')

mod='Cota20150604_1' # define model name
read=True ########### I N P U T #########################
readgeo=True ########### I N P U T #########################
geo_fname='grd.dat'
readdat=True ########### I N P U T #########################
dat_fname='flow2.inp'
readresults=False ########### I N P U T #########################
results_fname='flow2.out'
readflowH=True ########### I N P U T #########################
flowH_fname='FLOH.pkl'
readflowF=True ########### I N P U T #########################
flowF_fname='FLOF.pkl'

save=False ########### I N P U T #########################
savevtk=False ########### I N P U T #########################

print 'model=',mod
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+mod)    
if not os.path.exists('results'): 
    os.makedirs('results')   
    
if read:
    if readgeo is True: 
        print 'Reading geometry from '+ geo_fname
        geo=mulgrid(geo_fname)
        geo.radial=False # set flag that geometry is not yet radial.
    if readdat is True: 
        print 'Reading input data from '+ dat_fname
        dat=t2data(dat_fname) 
    if readresults is True:
        print 'Reading results from '+ results_fname
        results=t2listing(results_fname)
    if readflowH is True and os.path.isfile(flowH_fname):
        print 'Reading saturation data from '+ flowH_fname
        flowH=ptg.load_obj(flowH_fname)
        times=ptg.load_obj('time.pkl')
    else: 
        print('CANT READ flow FILE......Continuing with flowH={}')
        flowH=times={}
    if readflowF is True and os.path.isfile(flowF_fname):
        print 'Reading saturation data from '+ flowF_fname
        flowF=ptg.load_obj(flowF_fname)
        #times=ptg.load_obj('time.pkl')
    else: 
        print('CANT READ flow FILE......Continuing with flowF={}')
        flowF={}
t1=time.clock()        
print 'time to read=',(t1-t0)  
      
ipt.icepost(mod, geom_data=geo, results=results, times=times, save=save, savevtk=savevtk, flows={'FLOH':flowH,'FLOF':flowF})
print 'time to run =', time.clock()-tinit