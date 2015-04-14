# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/briochh/.spyder2/.temp.py
"""

import numpy as np
import scipy.interpolate
import os
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

plt.close('all')

# change directory  
os.chdir(r'C:\Users\glbjch\Google Drive\Molly')
data=np.genfromtxt('1dmodeltimes.txt',dtype=None,skiprows=4,names='k,t,Q')

## read in data
#KTQ=np.loadtxt('Report.txt',delimiter=' ',skiprows=4, usecols=(3,4,5)) 
#
## reorganise into three separate arrays. 
#K,T,Q=KTQ[:,0],KTQ[:,1],KTQ[:,2]
K,T,Q=data['k'],data['t'],data['Q']
#
## define a grid based on conductivity and heat flux
X,Y=np.meshgrid(K,Q)
#
## grid the data onto the new gird
Ts = ml.griddata(K,Q,T,X,Y,interp='linear')
#
## plot 
plt.imshow(Ts,extent=[X.min(),X.max(),Y.min(),Y.max()],aspect='auto',origin='lower')
plt.colorbar()
plt.xlim((X.min(),X.max()))
plt.ylim((Y.min(),Y.max()))
plt.show()

plt.figure()
plt.pcolormesh(X,Y,Ts)
plt.colorbar()
plt.xlim((X.min(),X.max()))
plt.ylim((Y.min(),Y.max()))
plt.show()
