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

def unique_rows(data):
    """function to find the unique rows in an array"""
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])
    
t0=time.clock()
os.chdir("C:\Users\glbjch\Local Documents\Work\Modelling\Steffi_GRAV\dev_files")

#read in file
data=np.genfromtxt('WaterBalance.out',dtype=None,skiprows=2,usecols=(0,1,7),names='date,rain,rech')
dates=np.array([[int(dst) for dst in d.split('/')] for d in data['date']]) # split out date entries
rainrech=np.vstack((data['rain'],data['rech'])) # rain anf recharge in ML/day across whole island
daterainrech=np.hstack((dates,rainrech.T)) # complete date, rainfall and recharge array.
#minyear=dates[:,2].min() # first year
#maxyear=dates[:,2].max() # lastyear
#maxi=dates.shape[0]-1 
#maxj=maxyear-minyear+1

montotrain=0 # set accumulaters to zero
montotrech=0
mon,yr=dates[0,1],dates[0,2] # set start date

c=0 # set counter
area=200.0*200.0*2594 # for conversion into per m2 

mony=np.ascontiguousarray(dates[:,1:])
nummonths=unique_rows(mony).shape[0] # unique rows function need a contiguous array

montrechts=np.zeros((nummonths,5)) # set up array to contain monthly timeseries

##monthlysampling
for day in daterainrech: # for every day in the data set
    if day[1] == mon: # if we are in the same month add rain and recharge to acumulaters
        montotrain=montotrain+day[3]
        montotrech=montotrech+day[4]
    else: # if not then were are now in the next month and the previous totals need to be added to the monthly dataset
        # calculate secondsin month      
        d1=datetime.datetime(yr,mon,1,0,0,0) 
        d2=datetime.datetime(yr,mon,calendar.monthrange(yr,mon)[1],23,59,59)
        dt=(d2-d1).total_seconds()+1
        #add data to monthly timeseries
        montrechts[c]=[yr,mon,montotrain,montotrech,dt]
        #reset acumulater 
        montotrain,montotrech=day[3],day[4]
        mon=int(day[1]) # unpdate to new month and year
        yr=int(day[2])
        c=c+1 # update counter

## do the same for the final month        
d1=datetime.datetime(yr,mon,1,0,0,0)
d2=datetime.datetime(yr,mon,calendar.monthrange(yr,mon)[1],23,59,59)
dt=(d2-d1).total_seconds()+1       
montrechts[c]=[yr,mon,montotrain,montotrech,dt]

#convert all to Kg/s/m2
montrechts[:,3]=montrechts[:,3]*(1e6/area)/montrechts[:,4]
montrechts[:,2]=montrechts[:,2]*(1e6/area)/montrechts[:,4]
mura,mure=[np.mean(row)*(1e6/area/(24.*3600.)) for row in rainrech] #kg/s/m2
#normallise by divising by the means
norm_monthrech=np.vstack((np.divide(montrechts[:,3],mure),montrechts[:,4])).T
#save data.
hdr1='normalised recharge\tduration (s)\taverage recharge = '+ str(mure)+ ' kg/s/m2'
np.savetxt('norm_monthrech.txt',norm_monthrech,header=hdr1,delimiter='\t',fmt='%10.5g')

hdr2='year\tmonth\trainfall(kg/m2/s)\trecharge(kg/m2/s)\tseconds'
fmt='%4i\t%02i\t%e\t%e\t%e'
np.savetxt('montrechts.txt',montrechts,header=hdr2,fmt=fmt,delimiter='\t')

