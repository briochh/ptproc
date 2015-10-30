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
model='Cota20150810_1'
version=2
yrsec=365.25*24*3600
models=[model,model+'_ptb'+str(version)]#,model+'_rtn']
times={}
ts=np.zeros(1)
glaclim=[250.,2500.]
wd='C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'
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
    #if logt:
    #    tscale=np.log10(tscale) 
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
    qin=np.ma.masked_array(qts,[qts>0]) # mask where flow is positive (i.e. out of the model)
    qout=np.ma.masked_array(qts,[qts<0]) # mask where flow is negative (i.e. in to the model)
    
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
    #ssqtsout[ssqts<0]=0
    #ssqtsin[ssqts>0]=0
    
    if flow=='FLOH': # calculate meltrate etc.
        unit='W'
        meltratematrix=(qout.T/3.35E5) # kg/s/m2
        #meltratematrix[meltratematrix<0]=0 # kg/s/m2
        ssmelt=(ssqts/3.35E5)
        # change in meltrate
        deltameltrate=np.subtract(meltratematrix,ssmelt) # kg/s/m2 ~ mm/s
        #meltrate just within glacier            
        glacmeltrate=meltratematrix.T[(X>glaclim[0]) & (X<glaclim[1])].T  # kg/s/m2 ~ mm/s
        #change in meltrate within glacier 
        deltaglacmeltrate=deltameltrate.T[(X>glaclim[0]) & (X<glaclim[1])].T # kg/s/m2 ~ mm/s
        meltrate=np.zeros(len(tscale))
        glacArea=0
        i=0
        for t in tscale:
            for x,A,r in zip(X,Area,meltratematrix[i]):
                if r > 0 and x < glaclim[1] and x > glaclim[0]:
                    meltrate[i]= meltrate[i] + (r*A) # kg/s
                    glacArea=glacArea+A
            i=i+1  
        meltrate_mmpyr= (meltrate*yrsec)/((np.pi*(glaclim[1]**2))-(np.pi*(glaclim[0]**2))) # kg/yr/m2 ~ mm/yr 
    else: unit ='kg/s'
    
    ## a quick plot of flows into atmosphere at X and time.
    plt.figure()
    plt.pcolormesh(X,tscale,qts.T, rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.1e')
    cbar.set_label(flow + r' out of the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Flow ('+flow+') out of the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig(stitchpath+model+'_'+flow+'.pdf',dpi=400)
    
    plt.figure()
    plt.pcolormesh(X,tscale,np.subtract(qout.T,ssqts), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.1e')
    cbar.set_label('Change in '+ flow + r' out of the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Change in flow ('+flow+') out of the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig(stitchpath+model+'_'+'delout_'+flow+'.pdf',dpi=400)

    
    plt.figure()
    plt.pcolormesh(X,tscale,-np.subtract(qin.T,ssqts), rasterized=True,cmap='rainbow') #W or (ks/s) /m2
    cbar=plt.colorbar(format='%.1e')
    cbar.set_label('Change in '+ flow + r' in to the model ('+ unit + r'/m$^{2}$)')
    #plt.xlim(0,2500)
    plt.ylim(tscale.min(),tscale.max())
    plt.title('Change in flow ('+flow+') in to the model')
    plt.xlabel('Distance from axial centre (m)')
    plt.ylabel('Time (yrs)')
    if save:
        plt.savefig(stitchpath+model+'_'+'delin_'+flow+'.pdf',dpi=400)
    
    if save:
        if os.path.isfile(stitchpath+flow+'.pkl'):
            print(flow+' flow alreay pickled')
        else:
            ptg.save_obj(stitch,stitchpath+flow+'.pkl')  
        if os.path.isfile(stitchpath+'alltime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,stitchpath+'alltime.pkl')
        if os.path.isfile(stitchpath+'ptbtime.pkl'):
            print('time alreay pickled')
        else:
            ptg.save_obj(ts,stitchpath+'ptbtime.pkl')

plt.figure()
plt.pcolormesh(X,tscale,meltratematrix, rasterized=True,cmap='rainbow') # mm/s
cbar=plt.colorbar(format='%.1e')
cbar.set_label(r'Melt rate (kg/s/m$^{2}$)')
#plt.xlim((0,2500))
plt.ylim(tscale.min(),tscale.max())
plt.xlabel('Distance from centre axis (m)')
plt.ylabel('Time (yrs)')    
plt.title('Melt rate')
if save:
    plt.savefig(stitchpath+model+'_'+'meltrate_.pdf',dpi=400)  
   
plt.figure()
plt.pcolormesh(X,tscale,deltameltrate, rasterized=True,cmap='rainbow') # mm/s
cbar=plt.colorbar(format='%.1e')
cbar.set_label(r'Change in melt rate (kg/s/m$^{2}$)')
#plt.xlim((0,2500))
plt.ylim(tscale.min(),tscale.max())
plt.xlabel('Distance from centre axis (m)')
plt.ylabel('Time (yrs)') 
plt.title('Change in melt rate')
if save:
    plt.savefig(stitchpath+model+'_'+'meltrate_delta_.pdf',dpi=400)  

plt.figure()
plt.pcolormesh(X[(X>glaclim[0]) & (X<glaclim[1])],tscale,deltaglacmeltrate, rasterized=True,cmap='rainbow') # mm/s
cbar=plt.colorbar(format='%.1e')
cbar.set_label(r'Change in melt rate (kg/s/m$^{2}$)')
plt.xlim((0,2500))
plt.ylim(tscale.min(),tscale.max())
plt.title('Change in glacial melt rate')
if save:
    plt.savefig(stitchpath+model+'_'+'glacmeltrate_delta_.pdf',dpi=400) 
    
plt.figure()
plt.pcolormesh(X[(X>glaclim[0]) & (X<glaclim[1])],tscale,glacmeltrate, rasterized=True,cmap='rainbow', vmax=1.0e-3) # mm/s
cbar=plt.colorbar(format='%.1e')
cbar.set_label(r'Glacial melt rate (kg/s/m$^{2}$)')
#plt.xlim((0,250))
plt.ylim(tscale.min(),100)
plt.xlabel('Distance from centre axis (m)')
plt.ylabel('Time (yrs)')    
plt.title('Melt rate')
plt.tight_layout()
if save:
    plt.savefig(stitchpath+model+'_glac_meltrate_.pdf',dpi=400)  
    
plt.figure()
plt.plot(tscale,meltrate_mmpyr)
#plt.xlim(0, 30000)
plt.xlabel('Time (yrs)')
plt.ylabel('Average melt rate at glacial base (mm/yr)')
plt.title('Average basal meltrate')
if save:
    plt.savefig(stitchpath+model+'_'+'basalmelt.pdf')
    
plt.figure()
plt.semilogx(tscale,meltrate*yrsec)
plt.xlim(0.08, 100)
#plt.ylim(0,4.5e7) #!!!!!!!!
plt.xlabel('Time (yrs)')
plt.ylabel('Rate of mass loss from glacier (kg/yr)')
plt.title('rate of mass loss')
if save:
    plt.savefig(stitchpath+model+'_'+'melt_ts.pdf')
    
        
t1=time.clock()        
print 'time to read=',(t1-t0)  
