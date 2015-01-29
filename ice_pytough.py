# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 09:57:38 2015
\n
Module containing fuctions for building and running pytough icey models. \n

Reading results\n

@author: Brioch Hemmings 
"""
from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
from t2listing import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import bisect
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import scipy.constants
import shutil
import time
import os


def simple_readres( modelname, savevtk=False, tough2_input=None, geom_data=None, results=None):
    mod=modelname
    current_d=os.getcwd()
    print os.getcwd()
    os.chdir(mod)
    print os.getcwd()
    ## read output file
    t0=time.clock()
    if tough2_input is None:
       dat=t2data('flow2.inp') # tough2 input from input file
    else: dat=tough2_input

    if geom_data is None:
        geo=mulgrid('2dgrd.dat') # geometry from gempetry file
    else: geo=geom_data    

    if results is None:
        results=t2listing('flow2.out') # read output file
    
    grid=dat.grid # define input grid
    t1=time.clock()
    t=t1-t0
    print 'time2read .out=',t
    
    ## create directory for results
    if not os.path.exists('results'):
        os.makedirs('results')
    os.chdir('results')
    
    if savevtk:
       t0=time.clock()
       results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
       t1=time.clock()
       t=t1-t0
       print 'time2writevtks',t
    os.chdir(current_d)   
    return results