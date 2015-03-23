# -*- coding: utf-8 -*-
"""
Python script for comparing modeled gravity timeseries from 1D simulations
Created on Fri Jun 27 12:48:47 2014

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

plt.close('all')

save='yes'
ref_mod='20150304_1_rad_main'
mod='20150304_1_rad_main'
t0=time.clock()

num=5

os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/'+ref_mod+'/results')
ref_grav=np.vstack(np.loadtxt('axsym_int_microgal5.dat')).T
#ref_grav=np.array([[1.,5.,12.],[5.,8.,10.]])

os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/'+mod+'/results')
grav=np.vstack(np.loadtxt('axsym_int_microgal'+str(num)+'.dat')).T
#grav=np.array([[2.,4.,8.,15.],[8.,2.,8.,15.]])

if grav[0][-1] > ref_grav[0][-1]:
    ref_grav=np.concatenate((ref_grav,np.array([[grav[0][-1]],[grav[1][-1]]])),axis=1)
im1=plt.figure()
yrsec=3600*24*365.25
yrsec=1
#yrsec=1
plt.plot(ref_grav[0]/yrsec,ref_grav[1],'-',linewidth=2)
plt.plot(grav[0]/yrsec,grav[1],'-',linewidth=2)



times=np.array([])
gravdif=np.array([])
i=0
for t,g in zip(grav[0],grav[1]):
   #print 't='+str(t)
   #print 'g='+str(g)
   while i < len(ref_grav[0]):
      #print i
      if t >=ref_grav[0][i] and t <= ref_grav[0][i+1]:
         times=np.concatenate((times,[t]),axis=0)
         #print 'times='+str(times)
         grad=((ref_grav[1][i+1])-(ref_grav[1][i]))/((ref_grav[0][i+1])-(ref_grav[0][i]))
         # print 'grad='+str(grad)
         refgatt=ref_grav[1][i]+((t-ref_grav[0][i])*grad)
         #print 'refgatt='+str(refgatt)
         gravdif=np.concatenate((gravdif,[g-refgatt])
         ,axis=0)
         break
      else: i=i+1    
            
        
plt.plot(times/yrsec,gravdif,linewidth=2) 
plt.xlim((0,times.max()/yrsec))
plt.ylim((-60,100))
plt.ylabel(r'$\Delta g$ (microgal)',fontsize=18)
plt.xlabel('Time (years)',fontsize=18)  
plt.tick_params(axis='both',labelsize=18)
if save is 'yes':                 
   np.savetxt('gravdiff'+str(num)+'.dat',zip(times,gravdif))
   im1.savefig('vsref'+str(num)+'.pdf')   
   im1.savefig('vsref'+str(num)+'.png')
   f = open('relresultxt'+str(num)+'.txt','w')
   f.write('Model = '+mod+'\n'
               'Base = 220150304_1_rad_main5 \n'
               'relgrav max (microgal) =' +str(gravdif.max())+'\n'
               'relgrav min (microgal) =' +str(gravdif.min())+'\n'
               'Max amplidute (relgrav)='+str(gravdif.max()-gravdif.min())+'\n')
   f.close()                      
            



