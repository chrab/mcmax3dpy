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
import time

class DataMCMax3D(object):
  """
  Data container for the outputs of an MCMax3D model.
  """
  def __init__(self,modelDir=".",outDir=None,name=""):
    """
    Parameters
    ----------
    modelDir : string
      The path to the main model directory
    outDir : string
      The path to the output directory of the model (DEFAULT: modelDir/output)      
    name : string
      A name of the model (DEFAULT: "").  

    
    Attributes
    ----------   
    """
    self.modelDir=modelDir    
    """ string :
    The directory of the model.
    """
    self.outDir=outDir
    """ string :
    The output directory of the model. DEFAULT: `None`
    """
    self._dataDir=modelDir
        
    self.name=name
    """ string :
    An optional name for the model. DEFAULT: empty string
    """
    self.zones=None
    """ array_like(:class:`mcmax3dpy.read.Zone`) :
    The zones of the model. 
    see :class:`mcmax3dpy.read.Zone` for details.
    """    
    self.MCseds=None
    """ array_like(:class:`mcmax3dpy.read.SED`) :
    The Monte Carlo Spectral Energy Distribution(s) for the models (SED) 
    see :class:`mcmax3dpy.read.SED` for details.
    """ 
    self.RTseds=None
    """ array_like(:class:`mcmax3dpy.read.SED`) :
    The raytracing Spectral Energy Distribution(s) for the models (SED) 
    see :class:`mcmax3dpy.read.SED` for details.
    """ 
    self.psizes=None   
    """ array_like(float,ndim=1) :
    array with the global dust sizes.
    """  
    
    if outDir is not None:
      self._dataDir=self._dataDir+"/"+outDir
      
      
  def __str__(self):
    output = "Info MCMax3D model \n"
    output += "Name: " + self.name
    output += "\nOutputdirectory: " + self._dataDir
    output += "\nZones: " + str(len(self.zones))
    output += " MCseds: " + str(len(self.MCseds))
    output += "\n"
    return output   
   
      

class Zone(object):
  """
  Data structure for an MCMax3D Zone.
   
  Currently only works for a spherical grid.
  
  Attributes
  ----------
  """
  def __init__(self):
    self.fname=None
    self.fname_abun=None
    self.minrho=1.e-22 #g/cm^-3
    # FIXME: 
    self.rhodparticle=1.98996067e0 #g/cm^3 
    """ float :
    The density of a dust particle. 
    TODO: still harcoded!!!
    """   
    self.nr=None
    """ int :
    Number of radial grid points.
    """   
    self.nt=None
    """ int :
    Number of theta grid points.
    """   
    self.np=None
    """ int :
    Number of phi grid points (azimuthal).
    """   
    self.nsize=None
    self.comp=None # the dust particles

    self.r=None
    """ array_like(float,ndim=3) :
    The radial grid points.
    `UNIT:` au, `DIMS:` (np,nt,nr)
    """   
    self.theta=None
    """ array_like(float,ndim=3) :
    The theta grid points.
    `UNIT:` rad, `DIMS:` (np,nt,nr)
    """   
    self.phi=None
    """ array_like(float,ndim=3) :
    The theta grid points.
    `UNIT:` rad, `DIMS:` (np,nt,nr)
    """   
    self.r_grid=None
    self.theta_grid=None
    self.phi_grid=None
        # 
    self.x0=None
    """ float :
    The x zero point of this zones.
    `UNIT:` au
    """       
    self.y0=None
    """ float :
    The y zero point of this zones.
    `UNIT:` au
    """       
    self.z0=None
    """ float :
    The z zero point of this zones.
    `UNIT:` au
    """       
    self.x=None
    """ array_like(float,ndim=3) :
    The cartesian x coordinates.
    `UNIT:` au, `DIMS:` (np,nt,nr)
    """       
    self.y=None
    """ array_like(float,ndim=3) :
    The cartesian y coordinates.
    `UNIT:` au, `DIMS:` (np,nt,nr)
    """       
    self.z=None
    """ array_like(float,ndim=3) :
    The cartesian z coordinates.
    `UNIT:` au, `DIMS:` (np,nt,nr)
    """       
  
    self.rhod=None
    """ array_like(float,ndim=3) :
    The dust density.
    `UNIT:` |gcm^-3|, `DIMS:` (np,nt,nr)    
    """
    self.rhog=None
    """ array_like(float,ndim=3) :
    The gas density.
    `UNIT:` |gcm^-3|, `DIMS:` (np,nt,nr)    
    """    
    self.rhogVer=None # Vertical column density integrataded from the top to the midplane of the disk
    self.rhodVer=None # Vertical column density integrataded from the top to the midplane of the disk
    self.temp=None
    """ array_like(float,ndim=3) :
    The temperature.
    `UNIT:` K, `DIMS:` (np,nt,nr)    
    """    
    self.chi=None
    """ array_like(float,ndim=3) :
    The UV radiation field.
    `UNIT:` Drain field, `DIMS:` (np,nt,nr)    
    """      
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
    
    self.sd=None
    """ array_like(float,ndim=2) :
    The surfacedensity read from the additional file. Contains the radius and the 
    surfacedensity.
    `UNIT:` , `DIMS:` (nr,2)    
    """
  
  # TODO maybe make the read routine a method function (like in the prodimo scripts)
  def read(self,infile,psizes=None):

    print("INFO: Read fits input ...")
    self.fname=infile
    
    self.fname_abun=self.fname.replace(".fits.gz","_abun.fits")
    
    fitsMCMax3D=fits.open(infile)
    fitsMCMax3D.info()
    
    self.x0=float(fitsMCMax3D[0].header["X0"])
    self.y0=float(fitsMCMax3D[0].header["Y0"])
    self.z0=float(fitsMCMax3D[0].header["Z0"])
    
    self.nr=int(fitsMCMax3D[0].header["NR"])
    # this is just for plotting (so that is looks nicer)
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

    # read the density
    self.rhod=fitsMCMax3D[4].data
    self.temp=fitsMCMax3D[5].data
    self.comp=fitsMCMax3D[6].data    
    # calculate a dummy dust size distribution
    
    self.rhog=fitsMCMax3D[7].data
    self.chi=fitsMCMax3D[11].data    
    if len(fitsMCMax3D)>12:
      self.AVrad=fitsMCMax3D[12].data
      
      
    if os.path.isfile(self.fname_abun):
      self.species,self.abundances=read_abun_fits(self.fname_abun)
    
    if psizes is not None:
      self.calc_amean(psizes)

    
    t = time.process_time()    
    print("INFO: Calculate vertical column densities ...")
    self.rhogVer=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.integrate_vertical(self.rhog,self.rhogVer)
    self.rhodVer=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.integrate_vertical(self.rhod,self.rhodVer)
    print("TIME: ",time.process_time() - t)
    
    fitsMCMax3D.close()
    
    
    # check for the surfacedensity
    # FIXME: not very nice
    sdfname=self.fname.replace(".fits.gz",".dat")
    sdfname=sdfname.replace("Zone","surfacedens")
    print("INFO: read "+sdfname)
    self.sd=numpy.loadtxt(sdfname)
    

    
    
  def _set_cartesian_coord(self):
 
    self.x=self.r*numpy.sin(self.theta)*numpy.cos(self.phi)
    self.y=self.r*numpy.sin(self.theta)*numpy.sin(self.phi)
    self.z=self.r*numpy.cos(self.theta)
    
    # also aplly the translation 
    self.x=self.x+self.x0
    self.y=self.y+self.y0
    self.z=self.z+self.z0
    
    return
  
  def calc_amean(self,psizes):
    """
    Calculate the meand dust size at each point in the disk 
    
    """
    
    print("INFO: Calculate dust moments ...")
    self.amean=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a1mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a2mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    self.a3mom=numpy.zeros(shape=(self.np,self.nt,self.nr))
    
    doit=(self.rhog >self.minrho)
    
    psizes2=psizes[:]**2.0
    psizes3=psizes[:]**3.0
    onethird=1.0/3.0
    
    for ip in range(self.np):
      for it in range(self.nt):
        for ir in range(self.nr):

          if not doit[ip,it,ir]:
            self.amean[ip,it,ir]=0.1 
            continue
          
          nparts=self.comp[0,:,0,ip,it,ir]/(self.rhodparticle*psizes3[:])
          sum_nparts= numpy.sum(nparts)
          mom1=numpy.sum(nparts*psizes)
          mom2=numpy.sum(nparts*psizes2)
          mom3=numpy.sum(nparts*psizes3)

          self.amean[ip,it,ir]=(mom3/sum_nparts)**(onethird) # this is the amean output from prodimo
          self.a1mom[ip,it,ir]=(mom1/sum_nparts)
          self.a2mom[ip,it,ir]=(mom2/sum_nparts)
          self.a3mom[ip,it,ir]=(mom3/sum_nparts)
          
            
          #print(self.amean[ip,ir,it])
            
          #print(self.amean[ip,ir,it])

  def integrate_vertical(self,intfield,outfield,chem_species=False):
    """
    Integrates a quantity (intfield) in vertical direction from the top
    to the midplane of the disk for each field in the grid. 
    
    Currently the implementation is rather approximate
    """
    
    pih=math.pi/2.0
    tocm=(1.0*u.au).cgs.value
    # first find the midplane (close to the midplane)
    # FIXME: check what happends if nt is odd    
    imid = int(self.nt/2-1)  
      
    thetagrid=numpy.abs(pih-self.theta[0,:,0])
    tanthetagrid=numpy.tan(thetagrid)
    # over this indices the integration is done, where it is assumed
    # that the theta grid is symetrical
    thetaidx=range(0,imid+1)
    
    for ir in range(self.nr):    
      # take the closest theta to the midplane
      rmid=self.r[0,imid,ir]*math.cos(thetagrid[imid])
      
      zgrid=rmid*tanthetagrid
              
      itprev=None
      irprev=None
      # calculate the highest point (avoid theta=0 and pi)    
      # now step to the next point    
      for it in thetaidx:
        r=math.sqrt(rmid**2.0+zgrid[it]**2.0)
        irz=numpy.argmin(numpy.abs(r-self.r[0,it,:]))
        
        # do the integration
        if itprev is not None: 
          dz=(zgrid[itprev]-zgrid[it])*tocm

          # indixed for other theta half
          itoh=self.nt-(it+1)
          itohprev=self.nt-(itprev+1)
          
          if chem_species:
            outfield[:,it,irz,:]=outfield[:,itprev,irprev,:]+0.5*(intfield[:,itprev,irprev,:]+intfield[:,it,ir,:])*dz
            outfield[:,itoh,irz,:]=outfield[:,itohprev,irprev,:]+0.5*(intfield[:,itohprev,irprev,:]+intfield[:,itoh,ir,:])*dz
          else:
            outfield[:,it,irz]=outfield[:,itprev,irprev]+0.5*(intfield[:,itprev,irprev]+intfield[:,it,ir])*dz
            outfield[:,itoh,irz]=outfield[:,itohprev,irprev]+0.5*(intfield[:,itohprev,irprev]+intfield[:,itoh,ir])*dz
            
          
        itprev=it
        irprev=irz

class SED(object):
  """ 
  Data container for one SED output of MCMax  
  """
  def __init__(self):
    self.wl=None
    """ array_like(float,dim=1) :
    The wavelenghts [micron]
    """
    self.nu=None
    
    self.fluxJy=None
    """ array_like(float,dim=1) :
    The flux of the whole domain [Jy]
    """
    self.fluxJyZones=None
    """ array_like(float,dim=2) :
    The fluxes for each individual zone. 
    """
    self.fluxJyStars=None
    """ array_like(float,dim=2) :
    The fluxes for each individual Star
    """


def read_abun_fits(fname):
  """
  Reads a fits file the contains the chemical abundances. This file is is produced by the 
  MCMax3D - ProDiMo interface. For more information contact Christian RAb.  
  

  Parameters
  ----------
  fname : str 
    the file nme of the corresponding fits 

  Returns
  -------
  
  species : array_like(str,ndim=1) :
    a list of all the chemical species names. 
    
  abun : array_like(float,ndim=4)
    Four dimensional array with the abundances of all chemical species. 
  """

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


def read_MCSpec(directory):  
  """
  Tries to read all Monte Carlo SEDs in the give directory and returns them as a list. 

  Parameters
  ----------
  directory : str 
    the directory to search for the SEDs 

  Returns
  -------
  array_like(:class:`mcmax3dpy.read.SED`)
    a list of SEDs 
  """
  fseds=sorted(glob.glob(directory+"/MCSpec*.dat"))
  if fseds is None: 
    return None

  print("INFO: read MCSpecs ...")

  seds=list()  
  for fsed in fseds:
    data=numpy.loadtxt(fsed)
    sed=SED()
    sed.wl=data[:,0]
    sed.fluxJy=data[:,1]
    seds.append(sed)

  return seds


def read_RTSpec(directory,nzones=1):  
  """
  Tries to read all Ray Tracing SEDs in the give directory and returns them as a list. 

  Parameters
  ----------
  directory : str 
    the directory to search for the SEDs 

  Returns
  -------
  array_like(:class:`mcmax3dpy.read.SED`)
    a list of SEDs 
  """
  fseds=sorted(glob.glob(directory+"/RTSpec*.dat"))
  if fseds is None: 
    return None

  print("INFO: read RTSpecs ...")

  seds=list()  
  for fsed in fseds:
    data=numpy.loadtxt(fsed)
    # case of only on wavelenght point, make ndim=2 array anyway
    if data.ndim==1:
      data=numpy.array([data])
    print(data,data.shape,data.ndim)
    sed=SED()
    sed.wl=data[:,0]
    sed.fluxJy=data[:,1]
    sed.fluxJyZones=data[:,2:2+nzones]
    sed.fluxJyStars=data[:,2+nzones:]
    seds.append(sed)

  return seds


def read(modelDir=".",outDir=None,readParticles=False):
  """
  Trys to read all the possible output of MCMax3D and puts it into one data structure
  (see :class:`mcmax3dpy.read.DataMCMax3D`)
  """
  
  print("INFO: Try to read everything ...")
  
  data=DataMCMax3D(modelDir,outDir)

  if readParticles:
    t = time.process_time()
    data.psizes=read_particle_sizes(modelDir+"/Particles")
    print("TIME: ",time.process_time() - t)
    
  
  data.zones=read_zones(data._dataDir,psizes=data.psizes)
  data.MCseds=read_MCSpec(data._dataDir)
  data.RTseds=read_RTSpec(data._dataDir,len(data.zones))
  
  return data
  

def read_particle_sizes(directory):
  """
  Reads the particle files. This is need to e.g. calculate the mand dust sizes etc. 
  
  FIXME: this routine can currently only deal with one particle type and one temperatur.
  
  
  Returns
  -------
  array_like(float,dim=1)
    Returns a list with the different particle sizes.
  
  """
  
  print("INFO: Read particle sizes ...")

  fnames=glob.glob(directory+"/particle0001_*_0001*.fits.gz")

  if fnames is None or len(fnames)==0:
    print("WARN: Could not read any particles in directory "+directory)
    return None
  
  fnames.sort()
   
  psizes=list()
  for fname in fnames:
    fitsf=fits.open(fname)
      
    psizes.append(float(fitsf[0].header["A1"]))
  
  return numpy.array(psizes)
  
  
            
def read_zones(directory,psizes=None):
  """
  Tries to read all Zones fits files in the given directory and returns the 
  results in a list. 
  
  Parameters
  ----------
  directory : str 
    the directory to search for the Zone files. 
    
  Returns
  -------
  array_like(:class:`mcmax3dpy.read.Zone`)
    List with all Zone object. 
  """
  
  print("INFO: read_zones ...")
  
  #do some stuff
  print(directory+"/Zone*.fits.gz")
  sortedZones=sorted(glob.glob(directory+"/Zone*.fits.gz"))
  print(sortedZones)
  
  zones=list()
  for zoneFile in sortedZones:
    zone=Zone()
    zone.read(zoneFile,psizes=psizes)
    
    # also check for the abundances
    # move the whole stuff in a separate routine and pass with the zone as 
    # argument    
    if os.path.isfile(zone.fname_abun):
      print("INFO: Read abundances "+zone.fname_abun)
      species,abun=read_abun_fits(zone.fname_abun)
      zone.species=species
      zone.abundances=abun
      
      zone.chem_cd=numpy.zeros(shape=(zone.np,zone.nt,zone.nr,len(zone.species)))
      nd=numpy.zeros(shape=(zone.np,zone.nt,zone.nr,len(zone.species)))
      # calculate the vertical column densities
      print("INFO: Calculate vertical column densities species ..")

      # number densities (broadcast does not work here)
      for i in range(len(zone.species)):
        nd[:,:,:,i]=zone.rhog/2.2844156663114814e-24*zone.abundances[:,:,:,i]
      t = time.process_time()  
      zone.integrate_vertical(nd,zone.chem_cd,chem_species=True)
      print("TIME integrate: ",time.process_time() - t)
      
    zones.append(zone)
    
  return zones



            