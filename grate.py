# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 14:20:42 2014
\n
Module containing fuction for calculating and plotting simulated rates of gravity change. \n
@author: Brioch Hemmings 
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import bisect
from scipy.interpolate import interp1d

mpl.rcParams['xtick.labelsize']=14

mpl.rcParams['ytick.labelsize']=14

mpl.rcParams['axes.labelsize']=14

def grate( in_ts, winlen, save="yes", input_in="yrs" ):
    """ plot_grate( timeseries, window_length, save_option, input_in )\n
    Calculate and plot simulated rates of gravity change. 
    """
    plottimes={}
    windg={}
    yrsec=3600*24*365.25
    if input_in is 'yrs': in_ts[:,0]=in_ts[:,0]*yrsec
    for win in winlen:
        n=0
        plott=np.empty(len(in_ts))
        dy=np.empty(len(in_ts))
        yrdg=np.empty(len(in_ts)) 
        for t,g in in_ts:
            te=t+win*yrsec 
            i=bisect.bisect(in_ts[:,0], te)
            if i == len(in_ts):
                dy[n] = in_ts[i-1,1]+((te-in_ts[i-1,0])*(in_ts[i-1,1]-in_ts[i-2,1])/(in_ts[i-1,0]-in_ts[i-2,0]))
            else:
                dy[n]=in_ts[i-1,1]+((in_ts[i,1]-in_ts[i-1,1])*(te-in_ts[i-1,0])/(in_ts[i,0]-in_ts[i-1,0]))
            yrdg[n]=dy[n]-g
            plott[n]=(t+win*yrsec/2)
            n+=1
        plottimes['win_'+str(win)]=plott
        windg['win_'+str(win)]=yrdg
    
    #times=plott/yrsec
        
    #im1=plt.figure()
    #plt.plot(times,yrdg)
    #plt.ylabel(r'$\Delta g$ (microgal)/'+str(win)+'yr')
    #plt.xlabel('Time (years)')
    #plt.axis([0.0, 110,None,None])
    #
    #im1.savefig('mugal_per_'+str(win)+'yr.pdf')  
    
    x = in_ts[:,0]/yrsec
    y = in_ts[:,1]
    f = interp1d(x, y)
    #f2 = interp1d(x, y, kind='cubic')
    
    
    
    xnew = np.linspace(x.min(), x.max(), 100000)
    ynew=f(xnew)
    dgbydt=np.gradient(ynew,np.diff(xnew)[0])
    
    output=np.array([['peryear',dgbydt.min(),dgbydt.max()]])
    
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
    plt.show()
    
    
    if save is 'yes':
           im2.savefig('mugal_all_per_yr.pdf')
           np.savetxt('minmaxs_test.txt',output,fmt='%s') 