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




plt.close('all')

save=True
ptype='rel'


t0=time.clock()
    
os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/')
if ptype=='rel':
    res=np.loadtxt('rel_minmax.txt', skiprows=1)
    res2=np.loadtxt('abs_minmax.txt', skiprows=1)
else:
    res=np.loadtxt('abs_minmax.txt', skiprows=1)


x=np.arange(len(res[0]))
color_idx = np.linspace(0, 5,6)
xlab=['flat1','flat2','topo1','topo2','flat3','flat4','topo3','topo4']

f=plt.figure()
n=0
l=[]

if ptype=='rel':
    for n,r in zip(color_idx,[[res2[i],res2[i+1]] for i in range(0,len(res2)-1,2)]):
        x=x+0.14
        with plt.style.context(('ggplot')):
            ci=np.array(matplotlib.colors.hex2color(plt.rcParams['axes.color_cycle'][int(n)]))
        #G=(ci.max()+ci.min())/2. # lightness
        #G=np.mean(ci) # average
        G=0.21*ci[0]+0.72*ci[1]+0.07*ci[2] #luminosity
        #c=np.array([e+(e*1000./100.) for e in ci])
        #c[c>1]=1.
        with plt.style.context(('ggplot')):
            plt.vlines(x,r[1],r[0], linewidth=6, linestyle='-',color=str(G))
        
x=np.arange(len(res[0])) # reset x for plotting           
for n,r in zip(color_idx,[[res[i],res[i+1]] for i in range(0,len(res)-1,2)]):
    x=x+0.14
    with plt.style.context(('ggplot')):
        l=l+[plt.vlines(x,r[1],r[0], linewidth=6, linestyle='-',color=plt.rcParams['axes.color_cycle'][int(n)])]

x=np.arange(len(res[0])) # reset x for plotting   
for n,r in zip(color_idx,[res[i] for i in range(len(res)) if i % 2 ==0]):
    #l=l+plt.plot(x,np.append(r,r[-1]),drawstyle='steps-post', linewidth=3, linestyle='--')
    with plt.style.context(('ggplot')):
        plt.hlines(r,x,x+1, linewidth=1, linestyle='--',color=plt.rcParams['axes.color_cycle'][int(n)])
    #n+=1    
#
#    

for n,r in zip(color_idx,[res[i] for i in range(len(res)) if i % 2 ==1]):
##    plt.plot(x,np.append(r,r[-1]),drawstyle='steps-post', linewidth=3,)
    with plt.style.context(('ggplot')):
         plt.hlines(r,x,x+1, linewidth=1, linestyle='-',color=plt.rcParams['axes.color_cycle'][int(n)])
#     



barlocs=[0,2,4,6]
#plt.minorticks_on
with plt.style.context(('ggplot')):
    plt.grid(b=False)   
    plt.ylim([-200,150])
    plt.bar(barlocs,[plt.yticks()[0][-1]-plt.yticks()[0][0]]*len(barlocs),bottom=plt.yticks()[0][0],color='lightgrey', width=1, edgecolor = "none")
    plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5],xlab )
    plt.xlabel('Model')
    plt.ylabel(r'Minimum and maximum $\Delta g$ ($\mu$Gal)')
    
    lgd=f.legend(l, ('P1','P2','P3','P4','P5','P6'),ncol=6, bbox_to_anchor=[0.5, 0.95],loc='center')
if save:
    f.savefig(ptype+'_steps.pdf',bbox_inches='tight')
