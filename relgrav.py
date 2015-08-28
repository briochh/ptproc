# -*- coding: utf-8 -*-
"""
Python script for comparing modeled gravity timeseries from 1D simulations
Created on Fri Jun 27 12:48:47 2014

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
import pytoughgrav as ptg

plt.close('all')

save='yes'
parent='20150814_2'
ref_mod='20150814_2_var1'
num=6
ref_ts='axsym_int_microgal'+str(num)+'.dat'
mod=ref_mod
t0=time.clock()

tss=['axsym_int_microgal'+str(num)+'.dat' for num in np.arange(1,7)]


os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/'+parent)

ptg.relgrav(ref_mod,mod,ref_ts,tss)


