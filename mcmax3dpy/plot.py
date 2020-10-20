'''
Created on 15 Nov 2017

@author: rab
'''
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches
import matplotlib.colors as matcolors
import mcmax3dpy.image as mimage
import astropy.wcs as wcs
import astropy.units as u
import copy

from scipy import ndimage


def _initfig(ax=None,projection=None,**kwargs):
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
      if projection is not None:
        fig, ax = plt.subplots(1, 1,figsize=_sfigs(**kwargs),
                               subplot_kw=dict(projection=projection))
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
  old_settings = np.seterr(divide='ignore') 
  array = np.log10(array)
  np.seterr(**old_settings)  # reset to default  
  return array


def plot_cuts_zones(zones,fieldname,centerZoneIdx=None,
                    vlim=[None,None],vlabel=None,clevels=None,patches=None,rlim=None,ip=0,
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
      vmin=np.log10(np.min(field))
    else:
      vmin=math.log10(vmin)
      
    if vmax is None:
      vmax=np.log10(np.max(field))
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
    #x=np.append(x,[np.zeros(shape=(zone.nr))],axis=0)
    #x=np.insert(x,0,[np.zeros(shape=(zone.nr))],axis=0)
    z=zone.z[ip,:,:]
    #z=np.append(z,[z[zone.nt-1,:]],axis=0)
    #z=np.insert(z,0,[z[0,:]],axis=0)
    val=fieldlog[ip,:,:]
    #val=np.append(val,[val[zone.nt-1,:]],axis=0)
    #val=np.insert(val,0,[val[0,:]],axis=0)
    
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
    #x=np.append(x,[np.zeros(shape=(zone.nr))],axis=0)
    #x=np.insert(x,0,[np.zeros(shape=(zone.nr))],axis=0)
    z=zone.z[int(zone.np/2)+ip,:,:]
    #z=np.append(z,[z[zone.nt-1,:]],axis=0)
    #z=np.insert(z,0,[z[0,:]],axis=0)
    val=fieldlog[int(zone.np/2)+ip,:,:]    
    #val=np.append(val,[val[zone.nt-1,:]],axis=0)
    #val=np.insert(val,0,[val[0,:]],axis=0)

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
    x=np.append(x,[x[0,:]],axis=0)
    y=zone.y[:,int(zone.nt/2),:]
    y=np.append(y,[y[0,:]],axis=0)
    val=fieldlog[:,int(zone.nt/2),:]
    val=np.append(val,[val[0,:]],axis=0)
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

def plot_midplane_zones(zones,fieldname,centerZoneIdx=None,
                    vlim=[None,None],vlabel=None,clevels=None,patches=None,rlim=None,
                    patchesAzimuthal=None,species=None,plotGrid=False,**kwargs):
  """
  Plots the the xy (rphi) planes considering all zones.
     
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
  
  print("plot_midplane_zones: ",fieldname)
  #ccolors=["black","blue","red","0.8"]
  ccolors=["0.7","0.6","0.5","0.4"]  
  ccolors=["0.4","0.4","0.4"]
  
  fig, ax = plt.subplots(1, 1)
  figsize = fig.get_size_inches() 
  #figsize = figsize*0.9
  #figsize[0] = figsize[0] * 1.58
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
      vmin=np.log10(np.min(field))
    else:
      vmin=math.log10(vmin)
      
    if vmax is None:
      vmax=np.log10(np.max(field))
    else:
      vmax=math.log10(vmax)
    
#    phi1=int(zone.np/4)  
#    phi2=int(zone.np/4*3)  

    levels = ticker.MaxNLocator(nbins=39).tick_values(vmax, vmin)    
    ticks = ticker.MaxNLocator(nbins=6, prune="both").tick_values(vmin, vmax)
  
    """
    Midplane cut
    """  
    # plot through the midplane (or close through the midplane
    # FIXME: this is just a workaround to make it look nicer (to fill the circle)
    # so the last point is the same as the first point in phi 
    x=zone.x[:,int(zone.nt/2),:]
    x=np.append(x,[x[0,:]],axis=0)
    y=zone.y[:,int(zone.nt/2),:]
    y=np.append(y,[y[0,:]],axis=0)
    # take an average, because it is not exaclty at zero. 
    # FIXME: check if I really have used the correc index
    val=(fieldlog[:,int(zone.nt/2),:]+fieldlog[:,int(zone.nt/2+1),:])/2.0
    val=np.append(val,[val[0,:]],axis=0)
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
   
  CB=fig.colorbar(CS3, ax=ax,pad=0.005,ticks=ticks,
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

def plot_radial_zone(zone,field,ylabel,ip=0,ylim=[None,None]):
  '''
  Plot a quantity in the midplane as function of radius for a single zone. 
  
  Parameters
  ----------

  ip : array_like(int), or int
    index of phi coordinate to plot. if array_like the radial plot is done 
    for all given ips.
  
  
  '''

  fig, ax = plt.subplots(1, 1)
    
  
  if type(ip) in [list,tuple]:
    ips=ip
  else:
    ips=[ip]
  
  
  #print(int(zone.nt/2-1),zone.z[ip,int(zone.nt/2-1),:])
  #print(int(zone.nt/2),zone.z[ip,int(zone.nt/2),:])
  #print(int(zone.nt/2+1),zone.z[ip,int(zone.nt/2+1),:])      

#   print(int(zone.np/2),zone.y[int(zone.np/2),int(zone.nt/2),:])
#   print(0,zone.y[0,int(zone.nt/2),:])
#  
#   print(int(zone.np/2+1),zone.y[int(zone.np/2+1),int(zone.nt/2),:])
#   print(1,zone.y[1,int(zone.nt/2),:])
# 
#   print(int(zone.np/2-1),zone.y[int(zone.np/2-1),int(zone.nt/2),:])
#   print(-1,zone.y[-1,int(zone.nt/2),:])



  for idxp in ips:
    # there is no z=0, so take the average from teh two closes ones.
    # I checed this for an even nt, these are the two closes ones to z=0
    # plot also the opposite site, should work for even grid numbers  
    y1=(field[idxp,int(zone.nt/2),:]+field[idxp,int(zone.nt/2-1),:])/2.0
    y2=(field[int(zone.np/2+idxp),int(zone.nt/2),:]+field[int(zone.np/2+idxp),int(zone.nt/2-1),:])/2.0
    y=np.concatenate((np.flip(y2),y1))
    
    # what the "radius" so simply always take x at 0 
    x1=zone.x[0,int(zone.nt/2),:]  
    x2=-x1
    #x2=zone.x[int(zone.np/2+idxp),int(zone.nt/2),:]  
    x=np.concatenate((np.flip(x2),x1))    
   
    # the phi_grid gives different values from phi, (I gues phi_grid are the borders. 
    # so better use phi here to be consistent
    ax.plot(x,y,label=str((zone.phi[idxp,int(zone.nt/2),0]*u.rad).to(u.deg)))
      
  if ylim[0] is not None or ylim[1] is not None:
    ax.set_ylim(ylim)
  
  ax.set_ylabel(ylabel)
  ax.set_xlabel("r [au]")
  #ax.semilogx()
  ax.semilogy()
  ax.legend()

#   pdf.savefig(transparent=False)
#   plt.close(fig)    
  return fig

  

def plot_sed(model,MC=True,RT=True,RTIdx=0,full=True,zones=True,zonesIdx=None,
             stars=True,starsIdx=None,**kwargs):
  """
  Plots the SED of the model.
  
  Parameters
  ----------
  model : :class:`mcmax3dpy.read.DataMCMax3D`
    The MCMax3D model
  
  zonesIdx : array_like(int,ndim=1)
    Indices of the zones to plot (starting from zero). 
    
  RTIdx : int
    which RT spectrum should be plotted (useful if there are more than one.
    Default: 0
  
  """
  
  fig,ax=_initfig(**kwargs)
  
  if model.MCseds is not None and MC is True and full is True:
    sed=model.MCseds[0]
    ax.plot(sed.wl,sed.fluxJy,label="MC")

  if model.RTseds is not None and RT==True:
    sed=model.RTseds[RTIdx]

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
  
  
def plot_image(image,field="I",vlims=None,extend=None,projection="wcs",xlims=None,ylims=None,vlog=True,
               xlabel="RA",ylabel="Dec",cblabel=None,pa=None,powerNormGamma=None,
               showGrid=False,**kwargs):
  '''
  
  Plots a single image produce by MCMax3D.
  
  Parameters
  ----------
  
  image : fits hdu
    something that is similar to what comes aut from method:`image.prep_image` 
  
  coords : str
    `wcs` use the wcs coordinate system for plotting
    `wcsrelative` use the wcs coordinate system put relative to the center
    `pixel` use the pixelcoordinate system
    clas:`astropy.wcs.WCS` directly use this object
    
  field : str
    can be either I, Q, U or PI, Qphi,Uphi . Default is I
    If field == ALL a 3 times 2 plot grid is shown, 

  
  '''
  
  labelfontsize=None
  if isinstance(projection,wcs.WCS):
    proj=projection       
  elif projection=="wcs":
    # naxis=2 beause we deal only with RA and Dec ... do not care about the other axis
    proj=wcs.WCS(image.header,naxis=2)
    xlabel="RA"
    ylabel="Dec"
    # FIXME: is not very flexible
    labelfontsize=6
  elif projection=="wcsrelative":
    # naxis=2 beause we deal only with RA and Dec ... do not care about the other axis
    proj=mimage.linear_offset_coords(wcs.WCS(image.header,naxis=2))
    xlabel="rel. RA [arcsec]"
    ylabel="rel. Dec [arcsec]"

    if xlims is not None:     
      xlims[0]=proj.wcs.crpix[0]+xlims[0]/np.abs(proj.wcs.cdelt[0])
      xlims[1]=proj.wcs.crpix[0]+xlims[1]/np.abs(proj.wcs.cdelt[0])    

    if ylims is not None:     
      ylims[0]=proj.wcs.crpix[0]+ylims[0]/np.abs(proj.wcs.cdelt[0])
      ylims[1]=proj.wcs.crpix[0]+ylims[1]/np.abs(proj.wcs.cdelt[0])    

    
  else: 
    proj=None
    xlabel="pixel"
    ylabel="pixel"
    
  fields=["I","Q","U","PI","Qphi","Uphi"]    
    
  if field=="ALL":
    subplot_kw=None
    if projection is not None:
      subplot_kw=dict(projection=proj)
  
    fig, axes = plt.subplots(2, 3,figsize=scale_figs((2.7,2.1)),subplot_kw=subplot_kw) 
    
    # flatten the list
    axes = [item for sublist in axes for item in sublist]
 
    fieldsloop=fields
    single=False
  else:
    fig,ax=_initfig(projection=proj,**kwargs)
    axes=[ax]
    fieldsloop=[field]
    single=True
  
  naxis=image.header["NAXIS"]

  for field,ax in zip(fieldsloop,axes):
    if naxis==2:
      values=image.data
    elif naxis == 3:   
      index=fields.index(field)
      #print("Plotting: "+fields[index])
      
      if index>3:
        Qphi,Uphi=mimage.combine_Polarizations(image.data[1],image.data[2],0)
        if field == "Qphi":
          values=Qphi
        elif field == "Uphi":
          values=Uphi
        else:
          print("Don not understand field: ",field)
      else:
        values=image.data[index]
    
    else:
      print("ERROR: Cannot deal with this image NAXIS>3") 
 
 
    if pa is not None: 
      # print("Rotate the image by: ",pa)
      # they way how the inclination in Images.out is defined requries to rotate the image by -90 (counterclockwise)
      # and then bye the PA also counterclockwise
      values=ndimage.rotate(values,-(90.0+pa),reshape=False,order=1)
    
    if vlog:
      with np.errstate(divide='ignore',invalid='ignore'):
        values=np.log10(values)
      
    if vlims is None:
      if vlog:
        vmax=np.max(values)+np.log10(0.8)
        vmin=np.max(values)-5.0
        extend="both"
      else:
        vmax=np.nanmax(values)
        vmin=np.nanmin(values)
        extend="neither"
    else:
      if vlog:
        with np.errstate(divide='ignore'):       
          vmin=np.log10(vlims[0])
          vmax=np.log10(vlims[1])
        extend="both"
      else:
        vmin=vlims[0]
        vmax=vlims[1]
        extend="both"
        
    if powerNormGamma is not None:
       norm=matcolors.PowerNorm(gamma=powerNormGamma)
    else:
      norm=None
    
    im = ax.imshow(values,vmin=vmin,vmax=vmax,origin="lower",cmap="inferno",norm=norm)
     
    if xlims is not None:
      ax.set_xlim(xlims)
  
    if ylims is not None:
      ax.set_ylim(ylims)
  
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)   
  
    # FIXME: also not very flexible, but it is also easy enough to do outside of the routine  
    if showGrid:
      ax.grid(color='white', ls='dashed',lw=0.2)  
    
    if labelfontsize is not None:
      ax.tick_params(axis='both', which='major', labelsize=labelfontsize) 
  
    #ax.contour(np.log10(image[zoomto:npix-zoomto,zoomto:npix-zoomto].value),levels=[0],colors="white",linewidths=1.0)
    cb=fig.colorbar(im,ax=ax,fraction=0.046, pad=0.01,extend=extend)
    if cblabel is None:
      if vlog:
        cblabel=r"log flux [$mJy\,arcsec^{-2}$]"
      else:
        cblabel=r"flux [$mJy\,arcsec^{-2}$]"
        
    cb.set_label(cblabel)
    ax.set_facecolor('black')
    for spine in ax.spines.values():
      spine.set_color('white')
  
    ax.tick_params(colors='white',labelcolor="black")
    
    if not single:
      # print the velocities relative to the systemic velocities
      props = dict(boxstyle='round', facecolor='white', edgecolor="none")
      ax.text(0.05, 0.95, field,transform=ax.transAxes, fontsize=6,fontweight="bold",
          verticalalignment='top', horizontalalignment="left", bbox=props)

  return fig  
