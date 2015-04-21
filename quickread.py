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


os.chdir('/Users/briochh/Documents/Workhere/testing')
mod='1Dtest2' # define model name
os.chdir(mod)
if not os.path.exists('results'): 
    os.makedirs('results')   
    
dat=t2data('flow2.inp')
grid=dat.grid
geo=mulgrid('grd.dat')
ptg.makeradial(geo,None,10.0)
results=t2listing('flow2.out')
os.chdir('results')
results.write_vtk(geo,'output.vtk',grid=grid,flows=True)
