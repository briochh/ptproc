# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 13:17:29 2014

Example file that calls and runs grate to calculate gravity rates for a given model simulation

@author: glbjch
"""

import os
os.chdir(r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough')
import pytoughgrav as ptg
import numpy as np
import matplotlib.pyplot as plt
import time

plt.close('all')
batch_or_straight='str' ########### I N P U T ######################### 
save=True ########### I N P U T #########################
nums=[1,2,3,4,5,6]

def anagrams(word):
    """ Generate all of the anagrams of a word. """ 
    if len(word) < 2:
        yield word
    else:
        for i, letter in enumerate(word):
            if not letter in word[:i]: #avoid duplicating earlier words
                for j in anagrams(word[:i]+word[i+1:]):
                    yield j+letter 

anas=[]
if __name__ == "__main__":
    for i in anagrams("batch"):
        anas=anas+[i]
        

if batch_or_straight in anas+['b','ba','bat','batc']:
    batch=True
    param='wt' ########### I N P U T #########################
    main=True ########### I N P U T #########################
else:
    batch=False
    mod='20150806_1_var' ########### I N P U T #########################

intype='' ########### I N P U T #########################
if intype is 'rel':
    infiles=['gravdiff'+str(n)+'.dat' for n in nums]#'axsym_int_microgal'+str(num)+'.dat' #'gravdiff1.dat')  ########### I N P U T #########################
else:
    infiles=['axsym_int_microgal'+str(n)+'.dat' for n in nums]#'axsym_int_microgal'+str(num)+'.dat' #'gravdiff1.dat')  ########### I N P U T #########################

input_times='yrs' #yrs ########### I N P U T #########################
windows=[2,5,10]  ########### I N P U T #########################

###########################################################################
############################# STRAIGHT MODE ###############################
if not batch:
    t0=time.clock()
    print 'running gravrates in straight mode (',batch_or_straight,')'
    print 'model=',mod
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Gravpaper/'+mod.split('_v')[0]+'/'+mod+'/results')
    #mod=mod+str(num)
    ptg.grate(mod,infiles,windows,input_in=input_times,save=save, intype=intype)



#################################################################
######################### BATCH MODE ############################
if batch:
    print 'running gravrates in batch mode (',batch_or_straight,')'
    print 'parameter=',param
    if main:
        print 'working on main model results'
    else:
        print 'working on warm up model results'
    os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
#input_in = 'yrs'
#save='yes'

#grav=vstack(np.loadtxt('microgal_5E-15.dat')).T
    if save:
        fall=open('gravrateslist.txt','w')
        fallmax=open('ratesmax.txt','w')
        fallmax.write('modelname \t ratemax1 \t ratemax2 \t ratemax5\t ratemax10  \n')
    for mod in [mod for mod in os.listdir('../'+param) if os.path.isdir(os.path.join('../'+param,mod))]:
        os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/batching/'+param)
        if main:
           os.chdir(mod+'/main/results')
        else:
           os.chdir(mod+'/results')
       # indata=np.loadtxt(infile)
        output=ptg.grate(mod,indata,windows,input_in=input_times,save=save,fall=fall,fallmax=fallmax)
    fall.close()
    fallmax.close()