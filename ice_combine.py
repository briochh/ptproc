# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:52:53 2016

@author: glbjch
"""

from t2data import *
from t2grids import *
from t2listing import *
import os
import pytoughgrav as ptg
import matplotlib.pyplot as plt
import numpy as np
import time
import ice_pytough as ipt
import copy
import matplotlib.colors
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec

class FixedOrderFormatter(mtick.ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        mtick.ScalarFormatter.__init__(self, useOffset=useOffset,
                                 useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag

t0=tinit=time.clock()
plt.close('all')

save=True ########### I N P U T #########################
#%% 
""" Define which models are to be combined in single plots.

Also define if models are 2D linear or axissymmetric 

and if the colour scales are loagarthmic """
modtype='FF'
model='Coto20150911_1_m2' # model - set log to false if need to
ptbrange=[1,2,3,4,5] # number of ptb models to plot
log=True # use params for log colour in LHS
#modtype='HPC1'
#model='Cota20150810_1_m2' # model - set log to false if need to
#ptbrange=[1,2,3,4] # number of ptb models to plot
#log=True # use params for log colour in LHS
#modtype='HPC2'
#model='Cota20150811_1_m2' # model - set log to false if need to
#ptbrange=[1,2,3,4] # number of ptb models to plot
#log=True # use params for log colour in LHS
#

#ptbcodes=['A','B','C','D','E']

yrsec=365.25*24*3600
models=['ptb'+str(ptb) for ptb in ptbrange]#,model+'_rtn'] # array of models to plot

times={}
ts=np.zeros(1)
glaclim=[0.,2500.] # define where the glacier is on the surface
wd='C:/Users/glbjch/Local Documents/Work/Modelling/Cotapaxi/'+model
os.chdir(wd)
flow='FLOH' # which flow to plot
times=ptg.load_obj('results/time.pkl') # load time data from inital model pickle
prelen=len(times) # length of inital model 

#%% set up figures
gs=gridspec.GridSpec(len(ptbrange), 2,wspace=0.02,width_ratios=[1,1]) # initalise grid for plots
combiplot=plt.figure(figsize=(8, 12)) 
fmatter=mtick.ScalarFormatter(useMathText=True, useOffset=False) # setup sci format for plots
fmatter.set_powerlimits((-2,4))


cbarax=combiplot.add_axes([0.13, 0.93, 0.35, 0.02]) # create color bar axis at top
combiplot200=plt.figure(figsize=(8, 12)) # zoomed in figure
cbarax200=combiplot200.add_axes([0.13, 0.93, 0.35, 0.02]) # create color bar axis at top
fmat='%.1'  # allow quick change of number format 

for ptb,ax in zip(models,ptbrange):
    mod=model+'_'+ptb
    stitch=ptg.load_obj('stitched_'+ptb+'/'+flow+'.pkl') # load stitched data
    ts=ptg.load_obj('stitched_'+ptb+'/alltime.pkl')
    if ts[0]==0.:
        ts=ts[1:]
    ptbtimes=ts[prelen:] # to just plot after perturbation
    tscale=np.subtract(ptbtimes,ptbtimes[0])/yrsec # convert to years
    X,Area,qts = ipt.prod_flowmatrix(stitch.values(),prelen) # generate x vs time matrix
    #craterflow=qts[X<glaclim[0]] # W/m2 or kg/s/m2
    glacflow=qts[(X<glaclim[1]) & (X>glaclim[0])]   
    #craterflowts,_ =ipt.convert2rate(qts.T,tscale,[0,250],Area,X)
    glacflowts,glac_contactArea =ipt.convert2rate(qts.T,tscale,glaclim,Area,X)
    if flow=='FLOH': # calculate meltrate etc.
        unit='W'
        #ipt.calcmeltrate(mod,qts,None,X,tscale,Area,
        #                 glaclim=glaclim,plot=False,save=save) 
    else: unit ='kg/s'
    
    
    #%%### LHS of full plot
    if log:
        vmax=None#4.0e3#1.0 # allow quick change colour limit
        vmin=1.0e-2
    else: 
        vmax=None#320
        vmin=0
        
    xlims=[0,2.5]
    tlims=[0,1000]
    lhax=combiplot.add_subplot(gs[(ax-1)*2],aspect=(max(xlims)-min(xlims))/(max(tlims)-min(tlims))) # create subplot for LHS

    if log:
        lhim=lhax.pcolormesh((X[(X<glaclim[1]) & (X>glaclim[0])])/1000.,
                tscale,np.ma.masked_array(glacflow,[glacflow<0]).T, 
                rasterized=True,cmap='rainbow',
                norm=matplotlib.colors.LogNorm(vmin=1e-2,vmax=vmax))
    else:
        lhim=lhax.pcolormesh((X[(X<glaclim[1]) & (X>glaclim[0])])/1000.,
                tscale,np.ma.masked_array(glacflow,[glacflow<0]).T, 
                rasterized=True,cmap='rainbow',vmin=0.0,vmax=vmax) # PLOT LHS
    lhax.set_xlim(xlims)
    lhax.set_ylim(tlims) # allow focus to particular timescale
    lhax.text(0.95,0.95,string.ascii_lowercase[(ax-1)*2]+')  '+modtype+'-'+string.ascii_uppercase[ax-1],verticalalignment='top', horizontalalignment='right',
        transform=lhax.transAxes,  style='italic')
    
    #%%RHS PLOT
    ax1=combiplot.add_subplot(gs[(ax*2)-1])
    plot_ax1, = ax1.plot(tscale,glacflowts)
    ax1.yaxis.set_major_formatter(fmatter)
    ax1.tick_params(axis='y', colors=plot_ax1.get_color())
    ax1.yaxis.label.set_color(plot_ax1.get_color())
    ax1.spines['left'].set_color(plot_ax1.get_color())
    ax1.yaxis.get_offset_text().set_color(plot_ax1.get_color())
    ax2=plt.twinx(ax1)
    plot_ax2, = ax2.plot(tscale,glacflowts/glac_contactArea,'g')
    #tlims=[0,1000]    
    
    if modtype=='FF':
        ylims2=[0, 250]
        ylims1=[0, 6.0e5]
    else:
        ylims2=[0,4.5]#[0, 3.5e-2]
        ylims1=[0,9.0e7]#[0, 7.0e5]
    #
    ax1.set_xlim(tlims)
    ax1.set_ylim(ylims1)
    ax2.set_ylim(ylims2)
    ax2.yaxis.set_major_formatter(fmatter)
    ax2.tick_params(axis='y', colors=plot_ax2.get_color())
    ax2.yaxis.label.set_color(plot_ax2.get_color())
    ax2.spines['right'].set_color(plot_ax2.get_color())  
    ax2.yaxis.get_offset_text().set_color(plot_ax2.get_color())
    ax2.text(0.95,0.95,string.ascii_lowercase[(ax*2)-1]+')  '+modtype+'-'+
            string.ascii_uppercase[ax-1],verticalalignment='top', 
            horizontalalignment='right',
            transform=ax2.transAxes,  style='italic')
        
    #%% ZOOM plot
    # LHS
    xlims=[0,2.5]
    tlims=[0,200]
    ax200=combiplot200.add_subplot(gs[(ax-1)*2],aspect=(max(xlims)-min(xlims))/(max(tlims)-min(tlims)))
    if log:
        im=ax200.pcolormesh(X[(X<glaclim[1]) & (X>glaclim[0])]/1000.,
                     tscale,np.ma.masked_array(glacflow,[glacflow<0]).T,
                     rasterized=True,cmap='rainbow',
                     norm=matplotlib.colors.LogNorm(vmin=1e-2,vmax=vmax)) # mm/s
    else:
        im=ax200.pcolormesh(X[(X<glaclim[1]) & (X>glaclim[0])]/1000.,
                     tscale,np.ma.masked_array(glacflow,[glacflow<0]).T,
                     rasterized=True,cmap='rainbow',vmin=vmin,vmax=vmax) # mm/s
    ax200.text(0.95,0.95, string.ascii_lowercase[(ax-1)*2]+')  '+modtype+'-'+string.ascii_uppercase[ax-1],verticalalignment='top', horizontalalignment='right',
        transform=ax200.transAxes,  style='italic')
    ax200.set_xlim(xlims)
    ax200.set_ylim(tlims)

    #%RHS
    if modtype=='FF':
        ylims2=[0, 250]
        ylims1=[0, 6.0e5]#7.0e5]
    else:
        ylims2=[0,4.5]#[0, 3.5e-2]
        ylims1=[0,9.0e7]#[0, 7.0e5]
    #
    ax2001=combiplot200.add_subplot(gs[(ax*2)-1])
    plot_ax1, = ax2001.plot(tscale,glacflowts)


    
    
#    #plt.ylim(craterflowts.min(),craterflowts.max()) #!!!!!!!!
#    
    ax2001.set_ylim(ylims1) 
    if ylims1[1] is not None:
        ax2001.yaxis.set_major_formatter(FixedOrderFormatter(np.floor(np.log10(ylims1[1])), useOffset=False, useMathText=True))
    else:
        ax2001.yaxis.set_major_formatter(fmatter)
    ax2001.yaxis.set_major_locator(mtick.MaxNLocator(nbins=5))
    #ax2001.yaxis.get_major_formatter(mtick.FormatStrFormatter('%0.0f'))
    ax2001.tick_params(axis='y', colors=plot_ax1.get_color())
    ax2001.yaxis.label.set_color(plot_ax1.get_color())
    ax2001.spines['left'].set_color(plot_ax1.get_color())
    ax2001.yaxis.get_offset_text().set_color(plot_ax1.get_color())
    ax2002=plt.twinx(ax2001)
    plot_ax2, = ax2002.plot(tscale,glacflowts/glac_contactArea,'g')
    ax2001.set_xlim(tlims)


    ax2002.set_ylim(ylims2)
    ax2002.yaxis.set_major_locator(mtick.MaxNLocator(nbins=5))
    ax2002.yaxis.set_major_formatter(fmatter)
    ax2002.tick_params(axis='y', colors=plot_ax2.get_color())
    ax2002.yaxis.label.set_color(plot_ax2.get_color())
    ax2002.spines['right'].set_color(plot_ax2.get_color())  
    ax2002.yaxis.get_offset_text().set_color(plot_ax2.get_color())
    ax2002.text(0.95,0.95,string.ascii_lowercase[(ax*2)-1]+')  '+modtype+'-'+string.ascii_uppercase[ax-1],verticalalignment='top', horizontalalignment='right',
        transform=ax2002.transAxes,  style='italic')
    # Dont want all axis labels
    if ax==ptbrange[-1] :  # if last row x labels
        lhax.set_xlabel('Distance x (km)')
        ax1.set_xlabel('Time (yrs)')
        ax200.set_xlabel('Distance x (km)')
        ax2001.set_xlabel('Time (yrs)')
    else:
        lhax.set_xticklabels([])
        #lhax.get_xaxis().set_visible(False)
        ax1.set_xticklabels([])
        ax200.set_xticklabels([])
        ax2001.set_xticklabels([])
        ax2002.set_xticklabels([])
    if len(ptbrange)%2==0 and ax/float(len(ptbrange))==0.5: #if an even number of ptbs then yaxis need to be offset between middle models    
 # if middle plot
        lhax.set_ylabel('Time (yrs)')
        lhax.yaxis.set_label_coords(-0.25, -0.1) 
        ax2.set_ylabel(r'Average heat-flux into glacier ('+ unit + r'm$^{-2}$)')
        ax2.yaxis.set_label_coords(1.15, -0.1) 
        ax1.set_ylabel(r'Total heat-flow into glacier ('+unit+')')
        ax1.yaxis.set_label_coords(-0.07, -0.1) 
        ax200.set_ylabel('Time (yrs)')
        ax200.yaxis.set_label_coords(-0.25, -0.1)
        ax2002.set_ylabel(r'Average heat-flux into glacier ('+ unit + r'm$^{-2}$)')
        ax2002.yaxis.set_label_coords(1.15, -0.1) 
        ax2001.set_ylabel(r'Total heat-flow into glacier ('+unit+')')
        ax2001.yaxis.set_label_coords(-0.07, -0.1)
    elif len(ptbrange)%2==1 and ax/float(len(ptbrange)+1)==0.5:
        lhax.set_ylabel('Time (yrs)')
        ax2.set_ylabel(r'Average heat-flux into glacier ('+ unit + r'm$^{-2}$)')
        ax1.set_ylabel(r'Total heat-flow into glacier ('+unit+')')
        ax200.set_ylabel('Time (yrs)')
        ax2002.set_ylabel(r'Average heat-flux into glacier ('+ unit + r'm$^{-2}$)')
        ax2001.set_ylabel(r'Total heat-flow into glacier ('+unit+')')
            


if log:
    cbar=combiplot.colorbar(lhim,cax=cbarax,orientation='horizontal',
                            ticks=mtick.LogLocator(subs=range(10)))
    cbar.set_label(r'Heat-flux into glacier ('+ unit + r'm$^{-2}$)', labelpad=-55)
    cbar200=combiplot200.colorbar(im, cax=cbarax200, orientation='horizontal',
                              ticks=mtick.LogLocator(subs=range(10)))#,format=fmat+'e')
    cbar200.set_label(r'Heat-flux into glacier ('+ unit + r'm$^{-2}$)', labelpad=-55)                                 
else:
    cbar=combiplot.colorbar(lhim,cax=cbarax,orientation='horizontal')
    cbar.set_label(r'Heat-flux into glacier ('+ unit + r'm$^{-2}$)', labelpad=-50)
    cbar.locator=mtick.MaxNLocator(nbins=5)
    cbar.update_ticks()
    cbar200=combiplot200.colorbar(im, cax=cbarax200, orientation='horizontal')
    cbar200.locator=mtick.MaxNLocator(nbins=5)
    cbar200.update_ticks()
    cbar200.set_label(r'Heat-flux into glacier ('+ unit + r'm$^{-2}$)', labelpad=-50) 


if save:
    combiplot.savefig(model+'combiplotplus.pdf',dpi=400,bbox_inches='tight')
    combiplot200.savefig(model+'combiplotplus200.pdf',dpi=400,bbox_inches='tight')

