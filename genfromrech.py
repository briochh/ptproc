# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 16:22:58 2015

@author: glbjch
"""

import time 
import shutil
import os
import numpy as np
import datetime
import calendar
import matplotlib.pyplot as plt

def unique_rows(data):
    """function to find the unique rows in an array"""
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])
    
t0=time.clock()
save=False
os.chdir("C:\Users\glbjch\Local Documents\Work\Modelling\Gravpaper\dev_files")

#read in file
data=np.genfromtxt('WaterBalance.out',dtype=None,skiprows=2,usecols=(0,1,7),names='date,rain,rech')
dates=np.array([[int(dst) for dst in d.split('/')] for d in data['date']]) # split out date entries
rainrech=np.vstack((data['rain'],data['rech'])) # rain anf recharge in ML/day across whole island
daterainrech=np.hstack((dates,rainrech.T)) # complete date, rainfall and recharge array.
#minyear=dates[:,2].min() # first year
#maxyear=dates[:,2].max() # lastyear
#maxi=dates.shape[0]-1 
#maxj=maxyear-minyear+1

t = np.array([datetime.datetime(date[2], date[1], date[0]) for date in dates])


montotrain=0 # set accumulaters to zero
montotrech=0
yrtotrain=0 # set accumulaters to zero
yrtotrech=0
mon,yr=dates[0,1],dates[0,2] # set start date

c=0 # set counter
y=0
area=200.0*200.0*2594 # for conversion into per m2 - 2594 is the number of land cells in recharge model 

mony=np.ascontiguousarray(dates[:,1:])
nummonths=unique_rows(mony).shape[0] # unique rows function need a contiguous array
numyrs=np.unique(mony[:,1]).shape[0] # unique rows function need a contiguous array
months=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
monthdic={month:[] for month in months}

montrechts=np.zeros((nummonths,5)) # set up array to contain monthly timeseries
yrrechts=np.zeros((numyrs,4)) # set up array to contain monthly timeseries
##monthlysampling
for day in daterainrech: # for every day in the data set
    if day[1] == mon: # if we are in the same month add rain and recharge to acumulaters
        montotrain=montotrain+day[3]
        montotrech=montotrech+day[4]     
        yrtotrain=yrtotrain+day[3]
        yrtotrech=yrtotrech+day[4]
    else: # if not then were are now in the next month and the previous totals need to be added to the monthly dataset
        # calculate secondsin month      
        d1=datetime.datetime(yr,mon,1,0,0,0) 
        d2=datetime.datetime(yr,mon,calendar.monthrange(yr,mon)[1],23,59,59)
        dt=(d2-d1).total_seconds()+1
        #add data to monthly timeseries
        montrechts[c]=[yr,mon,montotrain,montotrech,dt]
        monthdic[months[mon-1]].append(dt*montotrech/(24*3600))
        #reset acumulater 
        montotrain,montotrech=day[3],day[4]
        mon=int(day[1]) # unpdate to new month and year
        c=c+1 # update counter
        if day[2] == yr: 
            yrtotrain=yrtotrain+day[3]
            yrtotrech=yrtotrech+day[4]
        else: # if not then were are now in the next year and the previous totals need to be added to the yearly dataset
            # calculate secondsin yr     
            d1=datetime.datetime(yr,1,1,0,0,0) 
            d2=datetime.datetime(yr,12,31,23,59,59)
            dt=(d2-d1).total_seconds()+1
            #add data to yearly timeseries
            yrrechts[y]=[yr,yrtotrain,yrtotrech,dt]
            #reset acumulater 
            yrtotrain,yrtotrech=day[3],day[4]
            yr=int(day[2])
            y=y+1 # update counter




## do the same for the final month        
d1=datetime.datetime(yr,mon,1,0,0,0)
d2=datetime.datetime(yr,mon,calendar.monthrange(yr,mon)[1],23,59,59)
dt=(d2-d1).total_seconds()+1
dtyr=(d2-datetime.datetime(yr,1,1,0,0,0)).total_seconds()+1       
montrechts[c]=[yr,mon,montotrain,montotrech,dt]
monthdic[months[mon-1]].append(dt*montotrech/(24*3600))
yrrechts[y]=[yr,yrtotrain,yrtotrech,dtyr]

#convert all to Kg/s/m2
montrechts[:,3]=montrechts[:,3]*(1e6/area)/montrechts[:,4]
montrechts[:,2]=montrechts[:,2]*(1e6/area)/montrechts[:,4]
yrrechts[:,2]=yrrechts[:,2]*(1e6/area)/yrrechts[:,3]
yrrechts[:,1]=yrrechts[:,1]*(1e6/area)/yrrechts[:,3]
mura,mure=[np.mean(row)*(1e6/area/(24.*3600.)) for row in rainrech] #kg/s/m2
#normallise by divising by the means
norm_monthrech=np.vstack((np.divide(montrechts[:,3],mure),montrechts[:,4])).T
#save data.

tmon=np.array([datetime.datetime(np.int(date[0]),np.int(date[1]),1) for date in montrechts])
fig,ax1=plt.subplots()
ax1.bar(t,rainrech[1])
ax1.set_ylabel(r'Estimated whole island recharge (ML/day)')
ax1.set_xlabel('Time (years)')
ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
ax2=plt.twinx(ax1)
ax2.bar(tmon,montrechts[:,3],color='lightblue', edgecolor='lightblue', width=28)
ax2.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
ax2.set_ylabel(r'Estimated recharge rate (kg/s/m$^2$)')
ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
ax2.patch.set_visible(True) # hide the 'canvas' 
ax1.patch.set_visible(False) # hide the 'canvas' 

tyr=np.array([datetime.datetime(np.int(date[0]),1,1) for date in yrrechts])
fig2,ax1=plt.subplots()
ax1.bar(t,rainrech[1])
ax1.set_ylabel(r'Estimated whole island recharge (ML/day)')
ax1.set_xlabel('Time (years)')
ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
ax2=plt.twinx(ax1)
ax2.bar(tyr,yrrechts[:,2],color='lightblue', edgecolor='lightblue', width=365)
ax2.ticklabel_format(axis='y', style = 'sci', useOffset=False, scilimits=(-2,2))
ax2.set_ylabel(r'Estimated recharge rate (kg/s/m$^2$)')
ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
ax2.patch.set_visible(True) # hide the 'canvas' 
ax1.patch.set_visible(False) # hide the 'canvas' 

plt.figure()
monthstoplot=[]
med=[]
mean=[]
for mon in months:
    monthstoplot.append([monthdic[mon]])
    med.append([np.median(monthdic[mon])])
    mean.append([np.mean(monthdic[mon])])
boxes=plt.boxplot(monthstoplot,showmeans=True)
    

if save:    
    hdr1='normalised recharge\tduration (s)\taverage recharge = '+ str(mure)+ ' kg/s/m2'
    np.savetxt('norm_monthrech.txt',norm_monthrech,header=hdr1,delimiter='\t',fmt='%10.5g')

    hdr2='year\tmonth\trainfall(kg/m2/s)\trecharge(kg/m2/s)\tseconds'
    fmt='%4i\t%02i\t%e\t%e\t%e'
    np.savetxt('montrechts.txt',montrechts,header=hdr2,fmt=fmt,delimiter='\t')
    
    fig.savefig('Modeled_RechargeTS.pdf')

