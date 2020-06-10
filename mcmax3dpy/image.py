'''
Created on 5 Mar 2020

@author: rab
'''
import numpy as np
from astropy import wcs
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales



def combine_Polarizations(Q,U,phi0):
  '''
  Routine to calculate Q_phi and U_phi. 
  
  See for exampel Schmid+ 2006, Section 5 
  https://ui.adsabs.harvard.edu/abs/2006A%26A...452..657S/abstract
  '''
  
  Qphi= np.zeros(np.shape(Q))
  Uphi= np.zeros(np.shape(Q))
  x0= np.shape(Q)[1]/2#-0.5
  y0= np.shape(Q)[0]/2#-0.5
  phi0=(phi0*u.deg).to(u.rad).value
  for j in range(np.shape(Q)[0]): # over rows
    for i in range(np.shape(Q)[1]): # over columns
      phi= np.arctan2((float(i)-x0),(float(j)-y0))+phi0
      Qphi[i,j]= Q[i,j]*np.cos(2*phi)+U[i,j]*np.sin(2*phi)
      Uphi[i,j]= -Q[i,j]*np.sin(2*phi)+U[i,j]*np.cos(2*phi)
  
  return Qphi, Uphi


def prep_image(fname,fov,objectname=None,coordinates="00 42 30 +41 12 00",outfname=None,casaunits=False):
  '''
  Utility function to prepare a MCMax3D image fits file for e.g. plotting or for the use within CASA.
  
  Functionality 
  
  - adds a proper wcs coordinate system
  - adds the units for the fits file
  - can return or write the fits data to a fits file
  
  Parameters
  ----------
  fname : str 
    the path including the name to the fits file (can also by of type fits.gz)
    
  fov : float
    the fov of the original MCMax3D image in arcsec
    
  objectname : str
    a name for an astronomical object, the routine trys to finde the coordinates
    using :class:`astropy.coordinates.SkyCoord`
    
  coordinates : str
    corrdinsates as string in a format :class:`astropy.coordinates.SkyCoord` understands
    
  outfname : str
    if given the created fits file is written to a file with name given by `outfname`
    
  casaunits : boolean
    provide the fits file in units of JY/PIXEL (I think that is rather obsolete)
    
  Returns
  -------
    hdu : an `HDU` object class:`astropy.io.fits.PrimaryHDU`
  
  '''
  print("prep_image: "+fname)
  image=fits.open(fname)[0]  
  # first define the pixel scale
  npix=image.header["NAXIS1"]
  fov=fov*u.arcsec # arcsec
  pixelscale=(fov/npix)
  
  # assume that the center is really in the center
  centerpix=(npix-1)/2
  print("centerpix=",centerpix," pixelscale=",pixelscale)

  if objectname is not None:
    coord=SkyCoord.from_name(objectname)
  else:
    coord=SkyCoord(coordinates,unit=(u.hourangle, u.deg))    

  w = wcs.WCS(naxis=2)

  # Set up an "Airy's zenithal" projection
  # Vector properties may be set with Python lists, or Numpy arrays
  w.wcs.crpix = [centerpix, centerpix]
  # TODO: check if the minus is a general property or can depend on the actually object or image.
  # here I just used the example as given in the astropy documentaion 
  w.wcs.cdelt = np.array([-pixelscale.to(u.degree).value, pixelscale.to(u.degree).value])
  w.wcs.crval = [coord.ra.degree, coord.dec.degree]
  w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

  pixcrd = pixcrd = np.array([[1, 1], [centerpix, centerpix], [npix, npix]], dtype=np.float64)

  # Convert pixel coordinates to world coordinates.
  # The second argument is "origin" -- in this case we're declaring we
  # have 0-based (Numpy-like) coordinates.
  world = w.wcs_pix2world(pixcrd, 0)

  # Now, write out the WCS object as a FITS header
  header = w.to_header()

  header.append(("BTYPE","Intensity"))
  
  if casaunits:
    header.append(("BUNIT","JY/PIXEL"))
  else:
    header.append(("BUNIT","mJy/arcsec**2"))
  
  print(" ")
  #print(repr(header))

#   if pa > 0.0: 
#     for i in range(3):
#       print(i,image.data.shape)
#       image.data[i,:,:]=ndimage.rotate(image.data[i,:,:],pa,reshape=False)
  
  # header is an astropy.io.fits.Header object.  We can use it to create a new
  # PrimaryHDU and write it to a file.
  hdu = fits.PrimaryHDU(header=header)
  if casaunits:
    hdu.data=image.data[:,:,:]/1000.0*pixelscale.value**2      
  else:
    hdu.data=image.data[:,:,:]    
  
  if outfname is not None:
    hdu.writeto(outfname)

  return hdu


def linear_offset_coords(wcsabs, center=None):
  """
  Returns a locally linear offset coordinate system.
  
  Given a 2-d celestial WCS object and a central coordinate, return a WCS
  that describes an 'offset' coordinate system, assuming that the
  coordinates are locally linear (that is, the grid lines of this offset
  coordinate system are always aligned with the pixel coordinates, and
  distortions from spherical projections and distortion terms are not taken
  into account)
  
  taken from: https://github.com/aplpy/aplpy/issues/8
  
  Parameters
  ----------
  wcs : `~astropy.wcs.WCS`
      The original WCS, which should be a 2-d celestial WCS
  center : array_like
      The pixel coordinates on which the offset coordinate system should be
      centered. if `None` the centerpixels from wcsabs are taken
  """

  # Convert center to pixel coordinates
  #xp, yp = skycoord_to_pixel(center, wcs)
      
  # Set up new WCS
  
  new_wcs = wcs.WCS(naxis=2)
  #FIXME: do not know why +1 but otherwise the 0 points does not fit with the zero point in the pixel scale
  
  if center is None:
    new_wcs.wcs.crpix = wcsabs.wcs.crpix[0]+1, wcsabs.wcs.crpix[0]+1
  new_wcs.wcs.crval = 0., 0.
  # FIXME: assumes that the units are in degree
  scale=proj_plane_pixel_scales(wcsabs)*3600.
  # TODO: verify this somehow
  # I need to have the relative RA coordinates correct, I checked the absolute coordinats are projected correctly. 
  # This minus is also set in prep_image . I guess this has just to be consistent.
  scale[0]=scale[0]*-1.0 # otherwise it is the wrong direction, that does not change the orientation of the image just the labels
  new_wcs.wcs.cdelt = scale
  new_wcs.wcs.ctype = 'XOFFSET', 'YOFFSET'
  new_wcs.wcs.cunit = 'arcsec', 'arcsec'

  return new_wcs

