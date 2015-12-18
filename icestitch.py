# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 17:50:20 2015

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
import copy

#%% set up
t0=tinit=time.clock()
plt.close('all')
save=True ########### I N P U T #########################
model='Coto20150911_1_m2'
version=5
yrsec=365.25*24*3600
models=[model,model+'_ptb'+str(version)]#,model+'_rtn']
times={}
ts=np.zeros(1)
glaclim=[250.,2500.]
wd='C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'
#wd=os.path.expanduser("~")+'/GoogleDrive/Cotopaxi/'
#wd='C:/Users/glbjch/Local Documents/Work/Modelling/Molly project/'
if save:
    os.chdir(wd+model)
    if not os.path.exists('stitched_ptb'+str(version)): 
        os.makedirs('stitched_ptb'+str(version))
stitchpath=wd+model+'/stitched_ptb'+str(version)+'/'

flows=['FLOH','FLO(LIQ.)','FLO(GAS)']        
for flow in flows:
    os.chdir(wd)
    print flow
    data={}
    for mod in models:
        print 'model=',mod
        os.chdir(wd)
        os.chdir(mod)
        if ts[0]==0.0: # if first time through ts starts with 0
            times[str(mod)]=ptg.load_obj('results/time.pkl')
            ts=np.concatenate((ts,np.add(ts[-1],times[str(mod)])))
        data[str(mod)]=ptg.load_obj('results/'+flow+'.pkl')
    
    
    data['stitch']=copy.copy(data[str(models[0])])
    for conn in data['stitch']:
        for mod in models[1:]:
            q=data[str(mod)][conn][2]
            data['stitch'][conn]=data['stitch'][conn][0],data['stitch'][conn][1],np.concatenate((data['stitch'][conn][2],q))
    stitch=data['stitch']
    prelen=len(times[str(models[0])])
    if ts[0]==0.:
        ts=ts[1:]
    ptbtimes=ts[prelen:] # to just plot after perturbation
    tscale=np.subtract(ptbtimes,ptbtimes[0])/yrsec
    os.chdir(stitchpath)
    #if logt:
    #    tscale=np.log10(tscale) 
    # produce flow matrix (qts IN KG/S/M2 OR W/MS)
    X,Area,qts = ipt.prod_flowmatrix(stitch.values(),prelen) 
    #qin=qts.clip(max=0.)#np.ma.masked_array(qts,[qts>=0.]) # mask where flow is positive (i.e. out of the model)
    #qout=qts.clip(min=0.)#np.ma.masked_array(qts,[qts<=0.]) # mask where flow is negative (i.e. in to the model)

    X=[]
    ssqts=[]
    for x,a,q in data[str(models[0])].values():
        X.append(x)
        ssqts.append(q[-1])
    xinds=np.array(X).argsort()
    X=np.array(X)[xinds]
    ssqts=np.array(ssqts)[xinds] # J/s/m2 or kg/s/m2
    #ssqtsout=np.copy(ssqts)
    #ssqtsin=np.copy(ssqts)
    #ssqtsout=ssqts.clip(min=0)
    #ssqtsin=ssqts.clip(max=0)
    
    if flow=='FLOH': # calculate meltrate etc.
        unit='W'
        ipt.calcmeltrate(mod,qts,ssqts,X,tscale,Area,
                         glaclim=[250.,2500.],save=save) 
    else: unit ='kg/s'
    
    ## a quick plot of flows into atmosphere at X and time.
    ipt.plotflows(mod,flow,qts,X,tscale,Area,unit,ssqts,glaclim,save=save)
    
    if save:
        if os.path.isfile(flow+'.pkl'):
            print(flow+' flow alreay pickled')
        else:
            ptg.save_obj(stitch,flow+'.pkl')  
        if os.path.isfile(stitchpath+'alltime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,'alltime.pkl')
        if os.path.isfile('ptbtime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,'ptbtime.pkl')
    
        
t1=time.clock()        
print 'time to read=',(t1-t0)  
