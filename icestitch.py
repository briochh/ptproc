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
model='Cota20150611_1'
yrsec=365.25*24*3600
models=[model,model+'_ptb',model+'_rtn']
times={}
ts=np.zeros(1)
if save:
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+model)
    if not os.path.exists('stitched'): 
        os.makedirs('stitched')
flows=['FLOH','FLO(LIQ.)','FLO(GAS)']        
for flow in flows:
    print flow
    data={}
    for mod in models:
        print 'model=',mod
        os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+mod+'/results')
        if ts[0]==0.0: # if first time through ts starts with 0
            times[str(mod)]=ptg.load_obj('time.pkl')
            ts=np.concatenate((ts,np.add(ts[-1],times[str(mod)])))
        data[str(mod)]=ptg.load_obj(flow+'.pkl')
    
    
    data['stitch']=copy.copy(data[str(models[0])])
    for conn in data['stitch']:
        for mod in models[1:]:
            q=data[str(mod)][conn][2]
            data['stitch'][conn]=data['stitch'][conn][0],data['stitch'][conn][1],np.concatenate((data['stitch'][conn][2],q))
    stitch=data['stitch']
    prelen=len(times[str(models[0])])
    ptbtimes=ts[prelen:]
    
    X=[]
    qts=[]
    Area=[]
    for x,a,q in stitch.values():
        X.append(x)
        qts.append(q[prelen:])
        Area.append(a)
    inds=np.array(X).argsort()
    X=np.array(X)[inds]
    Area=np.array(Area)[inds]
    qts=np.array(qts)[inds] # J/s/m2 or kg/s/m2
    
    if flow=='FLOH': # calculate meltrate etc.
        unit='W'
    else: unit ='kg/s'
    if ts[0]==0.:
        ts=ts[1:]
    tscale=np.subtract(ptbtimes,ptbtimes[0])/yrsec
    #if logt:
    #    tscale=np.log10(tscale) 
    ## a quick plot of flows into atmosphere at X and time.
    plt.figure()
    plt.pcolormesh(X,tscale,qts.T, rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.0e')
    cbar.set_label(flow + r' out of the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Flow ('+flow+') out of the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig('../../'+model+'/stitched/'+flow+'.pdf',dpi=400)

    qout=np.copy(qts)
    qin=np.copy(qts)
    qout[qout<0]=0
    qin[qin>0]=0
    X=[]
    ssqts=[]
    for x,a,q in data[str(models[0])].values():
        X.append(x)
        ssqts.append(q[-1])
    inds=np.array(X).argsort()
    X=np.array(X)[inds]
    ssqts=np.array(ssqts)[inds] # J/s/m2 or kg/s/m2
    ssqtsout=np.copy(ssqts)
    ssqtsin=np.copy(ssqts)
    ssqtsout[ssqts<0]=0
    ssqtsin[ssqts>0]=0
    
    plt.figure()
    plt.pcolormesh(X,tscale,np.subtract(qout.T,ssqtsout), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.0e')
    cbar.set_label('Change in '+ flow + r' out of the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Change in flow ('+flow+') out of the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig('../../'+model+'/stitched/delout_'+flow+'.pdf',dpi=400)

    
    plt.figure()
    plt.pcolormesh(X,tscale,-np.subtract(qin.T,ssqtsin), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.0e')
    cbar.set_label('Change in '+ flow + r' in to the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Change in flow ('+flow+') in to the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig('../../'+model+'/stitched/delin_'+flow+'.pdf',dpi=400)
    
    if save:
        os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+model)
        if not os.path.exists('stitched'): 
            os.makedirs('stitched')   
        if os.path.isfile('stitched/'+flow+'.pkl'):
            print(flow+' flow alreay pickled')
        else:
            ptg.save_obj(stitch,'stitched/'+flow+'.pkl')  
        if os.path.isfile('stitched/alltime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,'stitched/alltime.pkl')
        if os.path.isfile('stitched/ptbtime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,'stitched/ptbtime.pkl')

        
t1=time.clock()        
print 'time to read=',(t1-t0)  
