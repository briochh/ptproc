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


os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Steffi_GRAV')
mod='20150327_1' # define model name
os.chdir(mod)
if not os.path.exists('results'): 
    os.makedirs('results')   
    
dat=t2data('flow2.inp')
grid=dat.grid
geo=mulgrid('grd.dat')
results=t2listing('flow2.out')
os.chdir('results')
results.write_vtk(geo,'output.vtk',grid=grid,flows=True)
