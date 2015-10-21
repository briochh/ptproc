# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 15:27:15 2015

@author: glbjch
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os
import time
import pytoughgrav as ptg
import numpy as np

plt.close('all')

save='no'
mod='20150814_4_var1'
num=6


t0=time.clock()
    
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/'+mod.split('_var')[0]+'/'+mod+'/results')
refmod=np.loadtxt('axsym_int_microgal'+str(num)+'.dat').T

f, axarr = plt.subplots(3, 2, sharex='col', sharey='row', figsize=[8,10])
f2, axarr2 = plt.subplots(3, 2, sharex='col', sharey='row', figsize=[8,10])

for station in xrange(1,num+1):
    print station
    bouggwt=np.loadtxt('bouguer_wt_'+str(station)+'.dat').T
    bougsat=np.loadtxt('bouguer_sat'+str(station)+'.dat').T
    axisym=np.loadtxt('axsym_int_microgal'+str(station)+'.dat').T
    relgrav=np.loadtxt('gravdiff'+str(station)+'.dat').T
    times=axisym[0]
    xpo=int(np.remainder((station+1),2))
    ypo=int(np.ceil(station/2.)-1)
    print xpo,ypo
    l1,=axarr[ypo, xpo].plot(times,bouggwt[1],'k--',linewidth=2,label='Bouguer water table')
    l2,=axarr[ypo, xpo].plot(times,axisym[1],'k-',linewidth=2,label='Distibuted gravity') 
    l3,=axarr2[ypo, xpo].plot(refmod[0],refmod[1],'--',color='0.2',linewidth=1,label='Reference Station')
    l4,=axarr2[ypo, xpo].plot(times,axisym[1],'-',color='0.2',linewidth=1,label='Benchmark P'+str(station)) 
    l5,=axarr2[ypo, xpo].plot(relgrav[0],relgrav[1],'k-',linewidth=2,label='Relative gravity')
    axarr[ypo, xpo].text(5,150,'P'+str(station), fontsize=14,weight='bold')
    axarr2[ypo, xpo].text(5,90,'P'+str(station), fontsize=14,weight='bold')
    #axarr[ypo, xpo].set_xlabel('Time (years)',fontsize=16)
    #if xpo == 0 and ypo ==1:
     #   axarr[ypo, xpo].set_ylabel(r'$\Delta g$ (microgal)',fontsize=18)
    #.plot(times,bouggwt[1],'--k',times,axisym[1],'-k',linewidth=2,)
    
    #plt.xlabel('Time (years)',fontsize=16)
    axarr[ypo, xpo].axis([0.0, 100,-180,180])
    axarr2[ypo, xpo].axis([0.0, 100,-165,110])
xlab=f.text(0.5, -0.02, 'Time (years)', ha='center',fontsize=18)
ylab=f.text(-0.04, 0.5, r'$\Delta g$ ($\mu$gal)', va='center', rotation='vertical',fontsize=18) 
xlab2=f2.text(0.5, -0.02, 'Time (years)', ha='center',fontsize=18)
ylab2=f2.text(-0.04, 0.5, r'$\Delta g$ ($\mu$gal)', va='center', rotation='vertical',fontsize=18) 
lgd=f.legend((l1,l2), ('Bouguer water table','Distibuted gravity'),ncol=2, bbox_to_anchor=[0.6, 1.04],loc='center')
lgd2=f2.legend((l3,l4,l5), ('Reference Station','Benchmark P',r'Relative $\Delta g$'),ncol=3, bbox_to_anchor=[0.15, 1.04,0.85,0.01],loc='center',mode="expand", borderaxespad=0.)  
f.tight_layout()
f2.tight_layout()
f2.show()

f.savefig('grav_combined_'+mod+'.pdf', bbox_extra_artists=(lgd,ylab,xlab), bbox_inches='tight')
f2.savefig('relgrav_combined_'+mod+'.pdf', bbox_extra_artists=(lgd2,ylab2,xlab2), bbox_inches='tight')




