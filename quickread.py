# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/briochh/.spyder2/.temp.py
"""

from t2data import *
from t2grids import *
from t2listing import *
import os
import pytoughgrav as ptg


os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Steffi_GRAV')
mod='20150429_1' # define model name
os.chdir(mod)
if not os.path.exists('results'): 
    os.makedirs('results')   
    
dat=t2data('flow2.inp')
grid=dat.grid
geo=mulgrid('grd.dat')
width=geo.bounds[1][1]-geo.bounds[0][1] #10.0
ptg.makeradial(geo,None,width)
results=t2listing('flow2.out')
os.chdir('results')
results.write_vtk(geo,'output.vtk',grid=grid,flows=True)
