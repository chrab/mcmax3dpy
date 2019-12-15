'''
Created on 15 Nov 2017

@author: rab
'''
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals


import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches



def _initfig(ax=None,**kwargs):
    '''
    Inits Figure and Axes object. 
    
    If an Axes object is passed, it is returned together with the Figure object.
    
    This is for a single plot (i.e. only one panel)
    
    Returns
    -------
    :class:`~matplotlib.figure.Figure`
    :class:`~matplotlib.axes.Axes` 
    '''

    if ax is not None:
      fig=ax.get_figure()

    else:
      fig, ax = plt.subplots(1, 1,figsize=_sfigs(**kwargs)) 

    return fig,ax
  
  
def _sfigs(**kwargs):
  '''
  Scale the figure size from matplotlibrc by the factors given in the 
  array sfigs (in kwargs) the first element is for the width the second for
  the heigth
  '''            
  if "sfigs" in kwargs:
    fac=kwargs["sfigs"]
    return scale_figs(fac)
  else:
    return scale_figs([1.0,1.0])
  
  
def scale_figs(scale):
  '''
  Scale the figure size from matplotlibrc by the factors given in the 
  array scale the first element is for the width the second for
  the heigth.
  '''            
  figsize=mpl.rcParams['figure.figsize']
  
  return (figsize[0]*scale[0],figsize[1]*scale[1])

def _dokwargs(ax,**kwargs):
  '''
  Handles the passed kwargs elements (assumes that defaults are already set)
  TODO: make this a general function .... 
  '''
  if "ylim" in kwargs: 
    ax.set_ylim(kwargs["ylim"])
    
  if "xlim" in kwargs: 
    ax.set_xlim(kwargs["xlim"])
            
  if "xlog" in kwargs:
    if kwargs["xlog"]: 
      ax.semilogx()
    else:
      ax.set_xscale("linear")
    
  if "ylog" in kwargs:
    if kwargs["ylog"]: 
      ax.semilogy()
    else:              
      ax.set_yscale("linear")
    
  if "xlabel" in kwargs:
    ax.set_xlabel(kwargs["xlabel"])  

  if "ylabel" in kwargs:
    ax.set_ylabel(kwargs["ylabel"])
  
#   if self.title != None and (not "notitle" in kwargs):
#     if self.title.strip() != "":
#       ax.set_title(self.title.strip())
    
  if "title" in kwargs:
    if  kwargs["title"] != None and kwargs["title"].strip() != "":
      ax.set_title(kwargs["title"].strip())
    else:
      ax.set_title("")  


def plog(array):      
  """
  Just a utilty function to avoid error-messages when taking the log of an arry
  """
  # ignore divide by zero in log10
  old_settings = numpy.seterr(divide='ignore') 
  array = numpy.log10(array)
  numpy.seterr(**old_settings)  # reset to default  
  return array


def plot_cuts_zones(zones,fieldname,centerZoneIdx=None,
                    vlim=[None,None],vlabel=None,clevels=None,patches=None,rlim=None,
                    patchesAzimuthal=None,patchesVertical=None,species=None,plotGrid=False,**kwargs):
  """
  Plots the xz (rtheta) and the xy (rphi) planes considering all zones.
   
  Currently the vertical cut  (rphi) is made a phi=0 and phi=pi. But this is done for all zones.
  
  TODO: Check if phi=0 is the same in all zones (relative to each other). 
  TODO: provide parameters for the cuts (e.g. not a phi=0).
  TODO: Currently value is always plotted on a logscale
  
  Parameters
  ----------
  zones : array_like(ndim=1)
    A list of the zones that should be plotted.
  
  fieldname : string
    The data (3D structure) the should be printed (e.g. "temp", "chi", "rhod"). For details see 
    :class:`mcmax3dpy.read.Zone`
        
  centerZoneIdx : int
    Center the figure to the center of the zone with centerZoneIdx. DEFAULT: `None`
  
  vlim : array_like(ndim=1)
    The range `vlim=[vmin,vmax]` for the values (also used for the colorbar). DEFAULT: `[None,None]` 
    
  rlim : float
    The maximum "radius" around the center that should be shown. However, the plot is of course squared. 
    But this is usefull to zoom in on the center.

  """
  
  print("plot_cuts_zones: ",fieldname)
  #ccolors=["black","blue","red","0.8"]
  ccolors=["0.7","0.6","0.5","0.4"]  
  ccolors=["0.4","0.4","0.4"]
  
  fig, axes = plt.subplots(1, 2,sharex=False,sharey=False)
  figsize = fig.get_size_inches() 
  figsize = figsize*0.9
  figsize[0] = figsize[0] * 1.58
  #figsize[1] =figsize[0]*0.61839012926
  fig.set_size_inches(figsize)
  
  
  # determine the relative center
  x0=y0=z0=0
  if centerZoneIdx is not None:
    x0=zones[centerZoneIdx].x0
    y0=zones[centerZoneIdx].y0
    z0=zones[centerZoneIdx].z0  
      
  for zone in zones:
    
    field=getattr(zone, fieldname)    
    if field is None: continue
    if species is not None:
      field=field[:,:,:,zone.species.index(species)]
    
    fieldlog=plog(field)
  
    vmin=vlim[0]
    vmax=vlim[1]
    
    if vmin is None:
      vmin=numpy.log10(numpy.min(field))
    else:
      vmin=math.log10(vmin)
      
    if vmax is None:
      vmax=numpy.log10(numpy.max(field))
    else:
      vmax=math.log10(vmax)
    
#    phi1=int(zone.np/4)  
#    phi2=int(zone.np/4*3)  

    levels = ticker.MaxNLocator(nbins=39).tick_values(vmax, vmin)    
    ticks = ticker.MaxNLocator(nbins=6, prune="both").tick_values(vmin, vmax)
  
  
    """
    Vertical cut
    """
  # get the coordinates
  # fix the coordinates for nicer plotting    
    x=zone.x[0,:,:]  
    #x=numpy.append(x,[numpy.zeros(shape=(zone.nr))],axis=0)
    #x=numpy.insert(x,0,[numpy.zeros(shape=(zone.nr))],axis=0)
    z=zone.z[0,:,:]
    #z=numpy.append(z,[z[zone.nt-1,:]],axis=0)
    #z=numpy.insert(z,0,[z[0,:]],axis=0)
    val=fieldlog[0,:,:]
    #val=numpy.append(val,[val[zone.nt-1,:]],axis=0)
    #val=numpy.insert(val,0,[val[0,:]],axis=0)
    
    ax=axes[0]
    CS = ax.contourf(x-x0, z-z0, val,levels=levels,extend="both",zorder=-20)    
    for c in CS.collections:
      c.set_edgecolor("face")         
    ax.set_rasterization_zorder(-19)
  
    if clevels is not None:
      ax.contour(CS, levels=clevels, colors=ccolors, linestyles="-",linewidths=0.5)
      
     
    if plotGrid:
      ax.scatter(x-x0,z-z0,s=0.2,color="0.5")     
  
    # plot also the other side,although is likely not exactly on the x-axis    
    x=zone.x[int(zone.np/2),:,:]
    #x=numpy.append(x,[numpy.zeros(shape=(zone.nr))],axis=0)
    #x=numpy.insert(x,0,[numpy.zeros(shape=(zone.nr))],axis=0)
    z=zone.z[int(zone.np/2),:,:]
    #z=numpy.append(z,[z[zone.nt-1,:]],axis=0)
    #z=numpy.insert(z,0,[z[0,:]],axis=0)
    val=fieldlog[int(zone.np/2),:,:]    
    #val=numpy.append(val,[val[zone.nt-1,:]],axis=0)
    #val=numpy.insert(val,0,[val[0,:]],axis=0)

    CS2 = ax.contourf(x-x0, z-z0, val,levels=levels,extend="both",zorder=-20)       
      # This is the fix for the white lines between contour levels
    for c in CS2.collections:
      c.set_edgecolor("face") 
    ax.set_rasterization_zorder(-19)    
    ax.set_aspect("equal")
    ax.set_xlabel("y [au]",labelpad=0)
    ax.set_ylabel("z [au]",labelpad=0)
    if clevels is not None:
      ax.contour(CS2, levels=clevels, colors=ccolors, linestyles="-",linewidths=0.5,zorder=0)
  
    if patches is not None:
      for patch in patches:
        # use a copy of the patch so that it can be reused
        ax.add_patch(copy.copy(patch),match_original=True)
        
    if patchesVertical is not None:
      for patch in patchesVertical:
        ax.add_patch(copy.copy(patch))
        
    if rlim is not None:
      ax.set_xlim([-rlim,rlim])
      ax.set_ylim([-rlim,rlim]) 
      #ax.set_aspect(2)    
        
    if plotGrid:
      ax.scatter(x-x0,z-z0,s=0.2,color="0.5")
        

    """
    Midplane cut
    """  
    # plot through the midplane (or close through the midplane
    # FIXME: this is just a workaround to make it look nicer (to fill the circle)
    # so the last point is the same as the first point in phi 
    x=zone.x[:,int(zone.nt/2),:]
    x=numpy.append(x,[x[0,:]],axis=0)
    y=zone.y[:,int(zone.nt/2),:]
    y=numpy.append(y,[y[0,:]],axis=0)
    val=fieldlog[:,int(zone.nt/2),:]
    val=numpy.append(val,[val[0,:]],axis=0)
    ax=axes[1]
    CS3 = ax.contourf(x-x0,y-y0, val,levels=levels,extend="both",zorder=-20)      
      # This is the fix for the white lines between contour levels
    for c in CS3.collections:
      c.set_edgecolor("face")
    ax.set_rasterization_zorder(-19)        
      
    if clevels is not None:
      ax.contour(CS3, levels=clevels, colors=ccolors, linestyles="-",linewidths=0.5,zorder=0)
      
    ax.set_aspect("equal")
    ax.set_xlabel("x [au]",labelpad=0)
    ax.set_ylabel("y [au]",labelpad=0)
    
    if patches is not None:
      for patch in patches:
        ax.add_patch(copy.copy(patch))

    if patchesAzimuthal is not None:
      for patch in patchesAzimuthal:
        ax.add_patch(copy.copy(patch))

  
    if rlim is not None:
      ax.set_xlim([-rlim,rlim])
      ax.set_ylim([-rlim,rlim])    
 
    if plotGrid:
      ax.scatter(x-x0,y-y0,s=0.2,color="0.5")
   
  plt.tight_layout()   
   
  CB=fig.colorbar(CS3, ax=axes.ravel().tolist(),pad=0.005,ticks=ticks,
                  format="%3.1f")
  if vlabel is not None:
    CB.set_label(vlabel)

  return fig


def plot_sd(zones,**kwargs):
  """
  Plots the surfacedensity of the Zones.
  
  """

  fig,ax=_initfig(**kwargs)


  for zone in zones:
    ax.plot(zone.sd[:,0],zone.sd[:,1],color="darkgray")

  _dokwargs(ax,**kwargs)
  
  return fig

  

def plot_sed(model,MC=True,RT=True,full=True,zones=True,zonesIdx=None,
             stars=True,starsIdx=None,**kwargs):
  """
  Plots the SED of the model.
  Currently the Monte Carlo SED is used. And simple the first one is taken.
  
  Parameters
  ----------
  model : :class:`mcmax3dpy.read.DataMCMax3D`
    The MCMax3D model
  
  zonesIdx : array_like(int,ndim=1)
    Indices of the zones to plot (starting from zero). 
    
  
  """
  
  fig,ax=_initfig(**kwargs)
  
  if model.MCseds is not None and MC is True and full is True:
    sed=model.MCseds[0]
    ax.plot(sed.wl,sed.fluxJy,label="MC")

  if model.RTseds is not None and RT==True:
    sed=model.RTseds[0]

    marker=None
    if len(sed.wl)==1: marker="+"

    if full is True:
      ax.plot(sed.wl,sed.fluxJy,label="RT",marker=marker)

    nzones=len(model.zones)
    if zones==True:
      for i in range(len(sed.fluxJyZones[0,:])):
        if zonesIdx is None or i in zonesIdx:
          ax.plot(sed.wl,sed.fluxJyZones[:,i],label="RT Zone "+str(i),
                  linestyle="--",marker=marker)
          
    if stars==True:
      for i in range(len(sed.fluxJyStars[0,:])):
        if starsIdx is None or i in starsIdx:
          ax.plot(sed.wl,sed.fluxJyStars[:,i],label="RT Star "+str(i),
                linestyle=":",marker=marker)
    
    
  ax.set_xlabel(r"wavelength [$\mu m$]")
  ax.set_ylabel(r"flux [Jy]")
  
  
  ax.semilogx()
  ax.semilogy()
  
  _dokwargs(ax,**kwargs)
  
  ax.legend()
  
  return fig
  
  
