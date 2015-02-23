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
import gzip


def simple_readres( modelname, savevtk=False, tough2_input=None, geom_data=None, results=None):
    """Easy reading of TOUGH2 results into pyTOUGH.  Option to write out vtk results files. Will return pyTOUGH results object if desired.
    
    Usage: 
    
    results=simple_readres(modelname, kwargs)
    
    modelname = directory name of model as string, e.g. "20150202_1"
    
    
    Optional kwargs:
    
    savevtk: flag to save vtk results files [default=False]
    
    tough2_input: provide the name of the tough2 input file [dafault="flow2.inp"]
    
    geom_data: provide name of pyTOUGH geometry data file ["grd.dat"]
    
    results: if results object already exists [None]
    
    """    
    mod=modelname # define modelname
    current_d=os.getcwd() # find current directory 
    print os.getcwd() 
    os.chdir(mod) # change to model directory
    print('Reading '+mod+' results from: ') #
    print os.getcwd()
    ## read output file
    t0=time.clock()
    if tough2_input is None:
       dat=t2data('flow2.inp') # tough2 input from input file
    else: dat=tough2_input

    if geom_data is None:
        geo=mulgrid('grd.dat') # geometry from geometry file
    else: geo=geom_data    

    if results is None:
        resfile='flow2.out'
#        resfile=gzip.open('flow2.out.gz','r').read()
        results=t2listing(resfile) # read output file
    else:
        results=results
    
    grid=dat.grid # define input grid
    t1=time.clock()
    t=t1-t0
    print 'time2read .out=',t
    
    ## create directory for results
    if not os.path.exists('results'):
        os.makedirs('results')
    os.chdir('results')
    
    ## save vtks...
    if savevtk:
       t0=time.clock()
       results.write_vtk(geo,mod+'_out.vtk',grid=grid,flows=True, time_unit='y')
       t1=time.clock()
       t=t1-t0
       print 'time2writevtks',t
    os.chdir(current_d)   
    return results