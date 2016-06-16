# -*- coding: utf-8 -*-
"""
Created on Fri May 20 13:03:16 2016

@author: glbjch
"""
import os 
import numpy as np

os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/'+
         'Gravpaper/20150806_3/20150806_3_var1')
         
gravts=np.loadtxt('results/axsym_int_microgal1.dat')
gravtimes=gravts[:,0]*365.25

#%% load recharge data
ts=np.loadtxt(r'C:\Users\glbjch\Local Documents\Work\Modelling\Gravpaper\dev_files\20150527_2_rand.dat')
## Define GENER block
#fpms=1.01308803322 # flux per meter squared
#%% evevation dependent params
fm=3.64742287695e-08
fc=2.15803271989e-06

yrsec=3600*24*365.25
elev=400.0
mult=0.7
length=100*yrsec
wavelength=1*yrsec
maxlength=3e5*yrsec

times=[0.0]+np.arange((wavelength/2),length,wavelength/2).tolist()+[length,maxlength]
numt=len(times)
mingen=2.0e-7

#%% run function
gxc=[]
for i in xrange(0,numt-3,2):
    highgx=((1+mult)*((elev*fm)+(fc)))*ts[i/2]
    if highgx < (1+mult)*mingen: highgx=(1+mult)*mingen
    lowgx=((1-mult)*((elev*fm)+(fc)))*ts[i/2]
    if lowgx < (1-mult)*mingen: lowgx=(1-mult)*mingen
    gxc=gxc+[lowgx,highgx] # recharge time series
    
times=np.divide(times[0:-2],3600.0*24.0)
gxc=np.array(gxc)
rech=np.multiply(gxc,24*3600)

t0=0.0
r0=rech[0]
t1=gravtimes[0]
acc=(t1-t0)*r0
func=(1-np.exp(-(t1-t0)/0.7))*np.exp(-(t1-t0)/40.0)
g=[2*np.pi*6.67e-11*acc*func]
print t0,t1,r0
print acc,"\n"
t0=t1
for t1 in gravtimes[1:]:
    acct=0.0
    print t0,t1
    ts=times[(times>=t0) & (times<=t1)]
    rs=rech[(times>=t0) & (times<=t1)]
    func=(1-np.exp(-(t1-t0)/0.7))*np.exp(-(t1-t0)/40.0)
    #acc=(ts[0]-t0)*r0
    #print t0,ts[0],r0,acc
    #t0=t1
    #print r0
    #r0=rs[0]
    for t,r in zip(ts[0:],rs[0:]):
        acct=acct+(t-t0)*r0
        print t0,t, r0,acct
        r0=r
        t0=t
   # r=gxc[inds[-1]]
    print t0,t1,r0
    acct=acct+(t1-t0)*r0
    g=g+[2*np.pi*6.67e-11*func*acct]
    print acct,"\n" 
    t0=t1    

   # r0=r 
    
    
#acs=0 # acumulations at each results time.....
