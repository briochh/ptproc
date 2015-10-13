# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:43:27 2015

@author: glbjch
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os
import time
import pytoughgrav as ptg
import numpy as np


plt.style.use('ggplot')

plt.close('all')

save=True
ptype='rel'


t0=time.clock()
    
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/')
res=np.loadtxt(ptype+'_minmax.txt', skiprows=1)
x=np.arange(len(res[0]))
color_idx = np.linspace(0, 5,6)
xlab=['fl1','fl2','tp1','tp2','fl3','fl4','tp3','tp4']

f=plt.figure()
n=0
l=[]

for n,r in zip(color_idx,[res[i] for i in range(len(res)) if i % 2 ==0]):
    #l=l+plt.plot(x,np.append(r,r[-1]),drawstyle='steps-post', linewidth=3, linestyle='--')
    plt.hlines(r,x,x+1, linewidth=1, linestyle='--',color=plt.rcParams['axes.color_cycle'][int(n)])
    #n+=1    
#
#    
#plt.gca().set_color_cycle(None)
for n,r in zip(color_idx,[res[i] for i in range(len(res)) if i % 2 ==1]):
##    plt.plot(x,np.append(r,r[-1]),drawstyle='steps-post', linewidth=3,)
     plt.hlines(r,x,x+1, linewidth=1, linestyle='-',color=plt.rcParams['axes.color_cycle'][int(n)])
#     
    
for n,r in zip(color_idx,[[res[i],res[i+1]] for i in range(0,len(res)-1,2)]):
    x=x+0.14
    l=l+[plt.vlines(x,r[1],r[0], linewidth=6, linestyle='-',color=plt.rcParams['axes.color_cycle'][int(n)])]


barlocs=[0,2,4,6]
#plt.minorticks_on
plt.grid(b=False)   
plt.ylim([-200,150])
plt.bar(barlocs,[plt.yticks()[0][-1]-plt.yticks()[0][0]]*len(barlocs),bottom=plt.yticks()[0][0],color='lightgrey', width=1, edgecolor = "none")
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5],xlab )
plt.xlabel('Model')
plt.ylabel(r'Minimum and maximum $\Delta g$ ($\mu$gal)')

lgd=f.legend(l, ('P1','P2','P3','P4','P5','P6'),ncol=6, bbox_to_anchor=[0.5, 0.95],loc='center')
if save:
    f.savefig(ptype+'_steps.pdf',bbox_inches='tight')
