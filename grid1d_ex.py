# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 10:29:19 2014

Script for executing grid1D function in batch mode

@author: glbjch
"""

import os
os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough')
import pytoughgrav as ptg
from t2data import *
import numpy as np

perm=5.0e-13
poro=0.34
rp={'type':11, 'parameters':[0.1,0.0,0.0,0.5,0.0,None,1.0]}
cp={'type':11, 'parameters':[0.0,-5000.0,0.001618,0.85,None,None,0.1]}
rechshift=0
wt=50.0
#poros=range(10,50,5)
perms=[1.0e-15,5.0e-15,1.0e-14,5.0e-14,1.0e-13,5.0e-13,1.0e-12,5.0e-12,1.0e-11]

# cprp alpha in MPA^-1
#      m
#      Slr
#      epsilon
cprps={'VS':[458.7,0.85,0.16,0.00239],
    'CHnv':[1.631,0.74,0.041,0.00288],
    'CHnz':[0.3140,0.38,0.11,0.150241],
     'PTn':[1.529,0.85,0.1,0.0288],
   'No5_2':[200.0,0.85,0.2,0.003179]}
   
#rechshifts=range(-14,12,2) #rechvariation in mm/day * 10

wts=range(20,160,20)
   
param='cprp'
for cprpname,cprplist in cprps.iteritems():#perm in perms:
#
   rp={'type':11, 'parameters':[cprplist[2],0.0,0.0,0.5,0.0,None,1.0]}
   cp={'type':11,'parameters':[0.0,-1.0*(1/(cprplist[0]*10**-6)),cprplist[3],cprplist[1],None,None,cprplist[2]]}
   
   os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/test')
   dat=t2data('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/base/flow2.inp')
   mod='1D_'+param+'_'+np.str(cprpname)
   #poro=poro/100.0
   if not os.path.exists(param+'/'+mod):
      os.makedirs(param+'/'+mod)

   ptg.grid1D(mod,dat,param,perm=perm,poro=poro,rp=rp,cp=cp,rechshift=rechshift,wt=wt)