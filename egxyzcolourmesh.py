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
os.chdir('/Users/briochh/GoogleDrive/Molly')
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
xi = np.linspace(K.min(), K.max(), np.unique(K).shape[0])
yi = np.linspace(Q.min(), Q.max(), np.unique(Q).shape[0])
#
## grid the data onto the new gird
Ts = ml.griddata(K,Q,T,X,Y,interp='linear')
Ti = ml.griddata(K,Q,T,xi,yi,interp='linear')
#
## plot 
plt.figure()
plt.xlim((X.min(),X.max()))
plt.ylim((Y.min(),Q.max()))
plt.imshow(Ti.data,origin='lower',aspect='auto',extent=[K.min(), K.max(),Q.min(), Q.max()])
plt.colorbar()
plt.show()

plt.figure()
plt.pcolormesh(xi,yi,Ti)
plt.colorbar()
plt.xlim((X.min(),X.max()))
plt.ylim((Y.min(),Y.max()))
plt.show()
