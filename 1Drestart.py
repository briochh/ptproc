# -*- coding: utf-8 -*-
"""1D restart"""
"""
Created on Sun Jun 08 15:36:21 2014

@author: glbjch
"""

from t2grids import *
from t2data import * # import classes and routines for creating TOUGH2 files
from t2incons import *
import matplotlib
import os
import random

os.chdir('C:/Users/glbjch/Local Documents/Work/Modelling/Pytough/')
mod='highhigh2'

basemod='highhigh'
if not os.path.exists(mod):
    os.makedirs(mod)

# read template file    
dat=t2data(basemod+'/flow2.inp')
geo=mulgrid(basemod+'/2dgrd.dat')
grid=dat.grid




# INCON
# change initial conditions from base model SAVE file
dat.incon.clear()
# Read from existing file
inc=t2incon(basemod+'/flow2.sav')
inc.write(mod+'/flow2.inc')

dat.parameter['max_timestep']=3.0e6
dat.parameter['print_interval']=3
dat.parameter['timestep']=[1000.0]
dat.output_times['time']=[1000.0,3600.0,8.6400e+04,3.1558e+07,3.1558e+08,3.1558e+09,3.1558e+10]
dat.output_times['num_times_specified']=7
dat.output_times['num_times']=7
dat.parameter['option'][12]=2

dat.generator.clear()

# Define GENER block
fpms=7.7354e-6 # flux per meter squared
fm=3.24e-8
fc=-7.199e-7
mingen=2.0e-7
cols=[col for col in geo.columnlist]
count=0

# time dependant generation
mult=0.65
yrsec=3600*24*365.25
sixmonth=yrsec/2
times=[0.0]+np.arange(sixmonth,(100*yrsec),sixmonth).tolist()+[(100*yrsec),1.0e13]
numt=len(times)
allgens=[]
dat.clear_generators()
gxc=[]
rand=[]
#for i in xrange(0,numt-3,2): rand=rand+[random.gauss(1,0.31)]:
rand= np.loadtxt(
r'C:\Users\glbjch\Local Documents\Work\Modelling\Pytough\20140613_2_py_it\rand.dat') # load random data file
# rand=((rand-mean(rand))/2)+mean(rand)
    
savetxt(mod+'/rand.dat',rand)
for col in cols:
    count=count+1
    lay=geo.column_surface_layer(col)
    blkname=geo.block_name(lay.name,col.name)
    gx=(grid.block[blkname].centre[2]*fm)+fc
    if gx < mingen: gx=mingen# for elevation dependant recharge!
    for i in xrange(0,numt-3,2):
       highgx=((1+mult)*((grid.block[blkname].centre[2]*fm)+(fc)))*rand[i/2]
       if highgx < (1+mult)*mingen: highgx=(1+mult)*mingen
       lowgx=((1-mult)*((grid.block[blkname].centre[2]*fm)+(fc)))*rand[i/2]
       if lowgx < (1-mult)*mingen: lowgx=(1-mult)*mingen
       gxc=gxc+[lowgx,highgx]
    #gxc=((numt-2)/2)*[lowgx,highgx]+[gx,gx]
    gxc=gxc+[gx,gx]   
    ex=numt*[1.0942e5]
    gxa=np.multiply(col.area,gxc).tolist()
    allgens.append(gxa)
    gen=t2generator(name=' q'+col.name,block=blkname,type='COM1',gx=None,ex=None,hg=None,fg=None, rate=gxa, enthalpy=ex, time=times,ltab=numt,itab=numt-1)
    #gen=t2generator(name=' q'+col.name,block=blkname,type='COM1', gx=gx*col.area, ex=1.0942e5)
    dat.add_generator(gen)
    dat.history_generator.append(blkname)
    
allgens=array(allgens)
gensum=sum(row[:] for row in allgens)
tforplot=[times[0]]
tforplot=np.append(tforplot,[2*[j] for j in times[1:-1]]+[times[-2]+yrsec])
tforplot=hstack(tforplot)
gforplot=[2*[j] for j in gensum[0:-1]]
gforplot=hstack(gforplot)
im1=plt.figure()
plt.plot(tforplot/yrsec,gforplot)
plt.ylabel('Generation rate ($kg/s/m^2$)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Time (years)')
plt.axis([0.0, tforplot.max()/yrsec,None,None])
savefig(mod+'/rech.pdf')
savetxt(mod+'/genertot.txt',gforplot)

geo.write(mod+'/2dgrd.dat') 
       
# write vtk of input information
grid.write_vtk(geo,mod+'/inparam.vtk',wells=True)
   
# write tough2 input file   
dat.write(mod+'/flow2.inp')



