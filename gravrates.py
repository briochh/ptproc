# -*- coding: utf-8 -*-
"""
Created on Wed Oct 08 15:13:25 2014
script to calculate and plot gravity change over a given window length.
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
import scipy.constants
import scipy.integrate as integrate
import bisect
from scipy.interpolate import interp1d

plt.close('all')
rcParams['xtick.labelsize']=14
rcParams['ytick.labelsize']=14
rcParams['axes.labelsize']=14

input_in = 'yrs'
save='yes'
os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\20140613_2_py_it\results')
#grav=vstack(np.loadtxt('microgal_5E-15.dat')).T
sixmonthts=np.loadtxt('gravdiff1.dat')
winlen=[2,5,10] #yrs
plottimes={}
windg={}
if input_in is 'yrs': sixmonthts[:,0]=sixmonthts[:,0]*yrsec
yrsec=3600*24*365.25
for win in winlen:
    n=0
    plott=np.empty(len(sixmonthts))
    dy=np.empty(len(sixmonthts))
    yrdg=np.empty(len(sixmonthts)) 
    for t,g in sixmonthts:
        te=t+win*yrsec 
        i=bisect.bisect(sixmonthts[:,0], te)
        if i == len(sixmonthts):
            dy[n] = sixmonthts[i-1,1]+((te-sixmonthts[i-1,0])*(sixmonthts[i-1,1]-sixmonthts[i-2,1])/(sixmonthts[i-1,0]-sixmonthts[i-2,0]))
        else:
            dy[n]=sixmonthts[i-1,1]+((sixmonthts[i,1]-sixmonthts[i-1,1])*(te-sixmonthts[i-1,0])/(sixmonthts[i,0]-sixmonthts[i-1,0]))
        yrdg[n]=dy[n]-g
        plott[n]=(t+win*yrsec/2)
        n+=1
    plottimes['win_'+str(win)]=plott
    windg['win_'+str(win)]=yrdg

times=plott/yrsec
    
#im1=plt.figure()
#plt.plot(times,yrdg)
#plt.ylabel(r'$\Delta g$ (microgal)/'+str(win)+'yr')
#plt.xlabel('Time (years)')
#plt.axis([0.0, 110,None,None])
#
#im1.savefig('mugal_per_'+str(win)+'yr.pdf')  

x = sixmonthts[:,0]/yrsec
y = sixmonthts[:,1]
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')



xnew = np.linspace(x.min(), x.max(), 100000)
ynew=f(xnew)
dgbydt=np.gradient(ynew,np.diff(xnew)[0])

output=array([['peryear',dgbydt.min(),dgbydt.max()]])

im2, axarr=plt.subplots(2,sharex=True)

data, =axarr[0].plot(x,y,'-',label='data')
#interp, = plt.plot(xnew,f(xnew),'b-',label='linear')
axarr[0].set_ylabel(r'$\Delta g$ ($\mu$gal)')
#plt.xlabel('Time (years)') 
#ax2=plt.twinx()
color_cycle=axarr[1]._get_lines.color_cycle
next(color_cycle)
grad, = axarr[1].plot(xnew,dgbydt,'-',label='gradient', color='0.6' )
#plt.legend(loc='best')
leghand=[data,grad]
leglab=['data','$\Delta g$/yr']
for win in winlen:
    temp=axarr[1].plot( plottimes['win_'+str(win)]/yrsec,windg['win_'+str(win)],'-', markersize=12,label=r'$\Delta g$/'+str(win)+'yrs')
    leghand=leghand+temp
    leglab=leglab+[r'$\Delta g$/'+str(win)+'yrs']
    output=np.concatenate((output,[[win,windg['win_'+str(win)].min(),windg['win_'+str(win)].max()]]))
axarr[1].set_ylabel(r'$\Delta g$ per time ($\mu$gal/time)')
plt.xlabel('Time (years)')    
data.axes.legend(leghand,leglab,loc='best',ncol=2,fontsize=12)
plt.axis([0.0, 100,None,None])
#plt.show()


if save is 'yes':  
       t0=time.clock()
       im2.savefig('mugal_all_per_yr.pdf')
       savetxt('minmaxs_test.txt',output,fmt='%s') 
#       f = open('minmaxs_test.txt','w')
#       f.write(output)
#       f.write('Model = '+mod+'\n'
#               'Mass in col max (kg) =' +str(well_water_mass.max())+'\n'
#               'Mass in col min (kg) =' +str(well_water_mass.min())+'\n'
#               'Max col amplidute (mass)='+str(well_water_mass.max()-well_water_mass.min())+'\n'
#               'grav (Boug) max (microgal) =' +str(microgal.max())+'\n'
#               'grav (Boug) min (microgal) =' +str(microgal.min())+'\n'
#               'Max (Boug) amplidute (grav)='+str(microgal.max()-microgal.min())+'\n'
#               'grav (int_axsym) max (microgal)='+str((gravt-gravt[0]).max())+'\n'
#               'grav (int_axsym) min (microgal)='+str((gravt-gravt[0]).min())+'\n'
#               'Max (int_axsym) amplitude (microgal)='+str((gravt-gravt[0]).max()-(gravt-gravt[0]).min())+'\n')
#       f.close()
#       savetxt('ztro.dat',zt_density_matrix)                
#       savetxt('waterweight'+str(wellno)+'.dat',zip(times,well_water_mass))
#       savetxt('microgal'+str(wellno)+'.dat',zip(times,microgal))
#       savetxt('axsym_int_microgal'+str(wellno)+'.dat',zip(times,(gravt-gravt[0])))
#       im.savefig('elemental_cont'+str(wellno)+'.pdf')       
#       im2.savefig('sl_t_profile'+str(wellno)+'test.png',dpi=300)
#       #im2.savefig('sl_t_profile'+str(wellno)+'.eps')
#       im1.savefig('microgal'+str(wellno)+'.pdf')
#       intgravplt.savefig('axsym_int_grav'+str(wellno)+'.pdf')
##
#     #savefig('sl_t_profile.pdf')
#       t1=time.clock()
#       t=t1-t0
#       print 'time2saveplot',t
