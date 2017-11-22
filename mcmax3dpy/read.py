'''
Created on 15 Nov 2017

@author: rab
'''
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

from astropy.io import fits
import astropy.units as u
import numpy
import glob
import os
import math

class Zone(object):
  """
  Data structure for an MCMax3D Zone.
   
  Currently only works for a spherical grid
  """
   
  def __init__(self):
    self.fname=None
    self.fname_abun=None
    self.minrho=1.e-22 #g/cm^-3
    self.rhodparticle=1.98996067e0 #g/cm^3    
    self.nr=None
    self.nt=None
    self.np=None
    self.nsize=None
    self.psizes=None # array of the mean dust radii for each size bin
    self.comp=None # the dust particles

    self.r=None
    self.theta=None
    self.phi=None

    self.r_grid=None
    self.theta_grid=None
    self.phi_grid=None
        # 
    self.x=None
    self.y=None
    self.z=None
    
    self.rhod=None
    self.rhog=None
    self.rhogVer=None # Vertical column density integrataded from the top to the midplane of the disk
    self.temp=None
    
    self.chi=None
    self.AVrad=None
    
    
    # a mean grain radius
    self.amean=None
    # dust size moments, used for the chemistry
    self.a1mom=None
    self.a2mom=None
    self.a3mom=None
    
    self.species=None
    self.abundances=None
    # the vertical column densities for the chemical species
    # shape(np,nt,nr,nspecies)
    self.chem_cd=None
  
  # TODO maybe make the read routine a method function (like in the prodimo scripts)
  def read(self,infile):

    print("INFO: Read fits input ...")
    self.fname=infile
    
    self.fname_abun=self.fname.replace(".fits.gz","_abun.fits")
    
    fitsMCMax3D=fits.open(infile)
    fitsMCMax3D.info()
    self.nr=int(fitsMCMax3D[0].header["NR"])
    # this is just for plotting (so that is looks nicer
    self.nt=int(fitsMCMax3D[0].header["NTHETA"])
    self.np=int(fitsMCMax3D[0].header["NPHI"])
    self.nsize=int(fitsMCMax3D[0].header["NSIZE"])

    self.r_grid=fitsMCMax3D[1].data    
    self.theta_grid=fitsMCMax3D[2].data
    self.phi_grid=fitsMCMax3D[3].data

    
    # hdu 0 is the whole grid in spherical coordinates
    # array ids 1. r,theta,phi 2. phi, 3. theta, 4. r
    self.r=fitsMCMax3D[0].data[0,:,:,:]
    self.r=fitsMCMax3D[0].data[0,:,:,:]
    # FIXME: work around to better match the outer radius in the plots
    # use the cell border for the outermost radius instead of the cell center
    self.r[:,:,-1]=((self.r_grid[-1]*u.cm).to(u.au)).value
    self.theta=fitsMCMax3D[0].data[1,:,:,:]
    self.phi=fitsMCMax3D[0].data[2,:,:,:]
    # convert to cartesian
    self._set_cartesian_coord()

        
    self.theta=fitsMCMax3D[0].data[1,:,:,:]
    self.phi=fitsMCMax3D[0].data[2,:,:,:]
    # convert to cartesian
    self._set_cartesian_coord()


    # read the density
    self.rhod=fitsMCMax3D[4].data
    self.temp=fitsMCMax3D[5].data
    self.comp=fitsMCMax3D[6].data    
    # calculate a dummy dust size distribution
    
    self.rhog=fitsMCMax3D[7].data
    self.chi=fitsMCMax3D[11].data    
    if len(fitsMCMax3D)>12:
      self.AVrad=fitsMCMax3D[12].data
    
    # FIXME: disable for the moment
    self.psizes=self.read_particle_sizes(self.nsize)
    self.calc_amean()
    
    print("INFO: Calculate vertical column densities ...")
    self.rhogVer=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.integrate_vertical(self.rhog,self.rhogVer)
    
    fitsMCMax3D.close()
    
  def read_particle_sizes(self,nsize):
    
    print("INFO: Read particles ...")
    psizes=numpy.zeros(shape=(nsize))
    # FIXME: fix fixed filenames
    dirn="Particles"
     
    for i in range(nsize):
      fname=dirn+"/particle0001_"+"{:04d}".format(i+1)+"_0001_f0.71_f0.29.fits.gz"
      try:        
        fitsf=fits.open(fname)
      except FileNotFoundError:
        fname="../"+dirn+"/particle0001_"+"{:04d}".format(i+1)+"_0001_f0.71_f0.29.fits.gz"
        fitsf=fits.open(fname)
        
      #fitsf.info()
      psizes[i]=float(fitsf[0].header["A1"])
    
    return psizes
    
  def _set_cartesian_coord(self):
 
    self.x=self.r*numpy.sin(self.theta)*numpy.cos(self.phi)
    self.y=self.r*numpy.sin(self.theta)*numpy.sin(self.phi)
    self.z=self.r*numpy.cos(self.theta)
    
    return
  
  def calc_amean(self):
    """
    Calculate the meand dust size at each point in the disk 
    
    Ist more or less a dummy, is not correct, now it is simply weighte by the number of particles
    
    FIXME: there are still hardcoded values
    """
    
    
    print("INFO: Calculate dust moments ...")
    self.amean=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a1mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a2mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a3mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    
    doit=(self.rhog >self.minrho)
    
    
    psizes2=self.psizes[:]**2.0
    psizes3=self.psizes[:]**3.0
    
    
    for ip in range(self.np):
      for it in range(self.nt):
        for ir in range(self.nr):
          sumnipart=0.0
          mom1=0.0
          mom2=0.0
          mom3=0.0
          if not doit[ip,it,ir]:
            self.amean[ip,it,ir]=0.1 
            continue          
          for ipart in range(self.nsize):
            # calculate the number of particles, assuming a fixed density 
            nipart=self.comp[0,ipart,0,ip,it,ir]/(self.rhodparticle*psizes3[ipart])
            sumnipart=sumnipart+nipart
            mom1=mom1+nipart*self.psizes[ipart]
            mom2=mom2+nipart*psizes2[ipart]
            mom3=mom3+nipart*psizes3[ipart]
          #print(sumnipart,ntotpart)     
          # FIXME: at the borders of the domain the dust density is very low so check for this
          # seams only to be a problem in the radial direction
#          if (self.rhog[ip,it,ir]/self.rhod[ip,it,ir])>1e17:
#            self.amean[ip,it,ir]=self.amean[ip,it,ir-1]
#            self.a1mom[ip,it,ir]=(mom1/sumnipart)
#            self.a2mom[ip,it,ir]=(mom2/sumnipart)
#            self.a3mom[ip,it,ir]=(mom3/sumnipart)
#          else:                 
          self.amean[ip,it,ir]=(mom3/sumnipart)**(1.0/3.0) # this is the amean output from prodimo
          self.a1mom[ip,it,ir]=(mom1/sumnipart)
          self.a2mom[ip,it,ir]=(mom2/sumnipart)
          self.a3mom[ip,it,ir]=(mom3/sumnipart)
          
            
          #print(self.amean[ip,ir,it])
            
          #print(self.amean[ip,ir,it])

  def integrate_vertical(self,intfield,outfield):
    """
    Integrates a quantity (intfield) in vertical direction from the top
    to the midplane of the disk for each field in the grid. 
    
    Currently the implementation is rather approximate
    """
    
    pih=math.pi/2.0
    # first find the midplane (close to the midplane)
    # FIXME: check what happends if nt is odd    
    imid = int(self.nt/2-1)  
      
    thetagrid=numpy.abs(pih-self.theta[0,:,0])
    # over this indices the integration is done, where it is assumed
    # that the theta grid is symetrical
    thetaidx=range(0,imid+1)
    
    for ir in range(self.nr):    
      # take the closest theta to the midplane
      rmid=self.r[0,imid,ir]*math.cos(thetagrid[imid])
      
      zgrid=rmid*numpy.tan(thetagrid)
              
      itprev=None
      irprev=None
      # calculate the highest point (avoid theta=0 and pi)    
      # now step to the next point    
      for it in thetaidx:
        r=math.sqrt(rmid**2.0+zgrid[it]**2.0)
        irz=numpy.argmin(numpy.abs(r-self.r[0,it,:]))
        
        # do the integration
        if itprev is not None: 
          dz=((zgrid[itprev]-zgrid[it])*u.au).cgs.value
          
          outfield[:,it,irz]=outfield[:,itprev,irprev]+0.5*(intfield[:,itprev,irprev]+intfield[:,it,ir])*dz
          # other theta grid half
          itoh=self.nt-(it+1)
          itohprev=self.nt-(itprev+1)
          outfield[:,itoh,irz]=outfield[:,itohprev,irprev]+0.5*(intfield[:,itohprev,irprev]+intfield[:,itoh,ir])*dz
          
        itprev=it
        irprev=irz


def read_abun_fits(fname):
  
  # primary hdu with the abundances
  hdulist = fits.open(fname)
  # table hdu with the species names
  
  hdulist.info()
  
  abun=hdulist[0].data
  # Species is a table data thing with only one column
  # needs to be a list at the moment
  species=list(hdulist[1].data.field(0))
  
  hdulist.close()
  
  return species,abun

            
def read_zones(directory):
  """
  Tries to read all Zones fits files in the given directory and returns the 
  results in a list. 
  """
  sortedZones=sorted(glob.glob(directory+"/Zone*.fits.gz"))
  print(sortedZones)
  
  zones=list()
  for zoneFile in sortedZones:
    zone=Zone()
    zone.read(zoneFile)
    
    # also check for the abundances
    # move the whole stuff in a separate routine and pass with the zone as 
    # argument    
    if os.path.isfile(zone.fname_abun):
      print("INFO: Read abundances "+zone.fname_abun)
      species,abun=read_abun_fits(zone.fname_abun)
      zone.species=species
      zone.abundances=abun
      
      zone.chem_cd=numpy.zeros(shape=(zone.np,zone.nt,zone.nr,len(zone.species)))
      # calculate the vertical column densities
      for ispec in range(len(species)):
        zone.integrate_vertical(zone.rhog/2.2844156663114814e-24*zone.abundances[:,:,:,ispec],
                       zone.chem_cd[:,:,:,ispec])
      
      
    zones.append(zone)
    
  return zones



            