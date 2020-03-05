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


def prep_image(fname,fov,objectname=None,coordinates="00 42 30 +41 12 00",outfname=None):
  '''
  Utility funtion to prepare a fits file for e.g. plotting.
  
  Functionoality 
  
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
    corrdinsates as string an format :class:`astropy.coordinates.SkyCoord` understands
    
  outfname : str
    if given the create fits is writen to a file with name given by `outfname`

  Returns
  -------
    hdu : an `HDU` object class:`astropy.io.fits.PrimaryHDU`
      or None if the fits file should be written to a file
  
  '''
  image=fits.open(fname)[0]  
  # first define the pixel scale
  npix=image.header["NAXIS1"]
  fov=fov*u.arcsec # arcsec
  pixelscale=(fov/npix)
  print(pixelscale)
  
  # assume that the center is really in the center
  centerpix=(npix-1)/2
  print(centerpix)

  if objectname is not None:
    coord=SkyCoord.from_name(objectname)
  else:
    coord=SkyCoord(coordinates,unit=(u.hourangle, u.deg))    

  w = wcs.WCS(naxis=2)

  # Set up an "Airy's zenithal" projection
  # Vector properties may be set with Python lists, or Numpy arrays
  w.wcs.crpix = [centerpix, centerpix]
  w.wcs.cdelt = np.array([pixelscale.to(u.degree).value, pixelscale.to(u.degree).value])
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
  header.append(("BUNIT","mJy/arcsec**2"))

  # header is an astropy.io.fits.Header object.  We can use it to create a new
  # PrimaryHDU and write it to a file.
  hdu = fits.PrimaryHDU(header=header)
  hdu.data=image.data[0,:,:]    
  
  if outfname is not None:
    hdu.writeto(outfilename)

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
  # FIXME: assumbes that the units are in degree
  scale=proj_plane_pixel_scales(wcsabs)*3600.
  #scale[0]=scale[0]*-1.0 # otherwise it is the wrong direction, that does not change the orientation of the image just the labels
  new_wcs.wcs.cdelt = scale
  new_wcs.wcs.ctype = 'XOFFSET', 'YOFFSET'
  new_wcs.wcs.cunit = 'arcsec', 'arcsec'

  return new_wcs

