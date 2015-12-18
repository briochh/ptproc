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

bmod=''#'Cota20150619_1/Cota20150619_1_ptb/'
mod='Coto20150911_1' # define model name
radial=False
eos=3
read=True ########### I N P U T #########################
readgeo=True ########### I N P U T #########################
geo_fname='grd.dat'
readdat=True ########### I N P U T #########################
dat_fname='flow2.inp'
readresults=True########### I N P U T #########################
results_fname='flow2.out'
readflow=True ########### I N P U T #########################
#flowH_fname='results/FLOH.pkl'
#readflowF=True ########### I N P U T #########################
#flowF_fname='results/FLOLIQ.pkl'

save=False ########### I N P U T #########################
savevtk=False ########### I N P U T #########################
if eos==1:
    flows={'FHEAT':{},'FLO(AQ.)':{},'FLO(GAS)':{}}
else: flows={'FLOH':{},'FLO(LIQ.)':{},'FLO(GAS)':{}}

print 'model=',mod
#os.chdir(os.path.expanduser("~")+'/GoogleDrive/Cotopaxi/'+bmod+'/'+mod)
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+bmod+'/'+mod)    
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
    if readflow is True:
        if os.path.isfile('results/time.pkl'):
            times=ptg.load_obj('results/time.pkl')
        else: 
        	print('CANT READ time FILE ......Continuing with times={}')
        	times={}
        for flow in flows.keys():
            flow_fname='results/'+flow+'.pkl'
            if os.path.isfile(flow_fname):
                print 'Reading saturation data from '+ flow_fname
                flows[flow]=ptg.load_obj(flow_fname)
            else: 
                print('CANT READ flow FILE ('+flow_fname+'......Continuing with flows['+flow+']={}')
                flows[flow]={}
    else: 
        print('Not reading '+flow+'.pkl.flows['+flow+']={}')
        flows[flow]={}
        

t1=time.clock()        
print 'time to read=',(t1-t0)  
      
ipt.icepost(mod, geom_data=geo,tough2_input=dat, radial=radial, results=results, times=times, save=save, savevtk=savevtk,flows=flows)
print 'time to run =', time.clock()-tinit