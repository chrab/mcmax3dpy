import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch


def read_opacities(rundir,Nzones,fvSI,fvC):

    ############################################################
    #
    # This rutine reads the opacities (extinction,absportion and
    # scattering) computed by MCMax3D. Inputs parameters are: 
    # "rundir" [str] path to the simulation's folder.
    # "Nzones" [int] number of zones defined in input.dat
    # "fvSI" [float] volume fraction of Silicates
    # "fvC" [float] volume fraction of Carbon 
    # The last two parameters are must be written as they appear
    # in the "particlexxx_xxx_xxx_fvSI_fvC.fits.gz" files.
    #
    # The output are three fits files with cantaining the 
    # mass-weighted exctinction, absorption and scattering 
    # opacities in cm^2/g
    #
    ############################################################


    ############################################################
    # Initializing matrices and arrays
    Nbins=np.zeros(Nzones)
    apows=np.zeros(Nzones)
    psizes_min=np.zeros(Nzones)
    psizes_max=np.zeros(Nzones)
    psizes=[]


    ############################################################
    # Reading information at rundir
    folder=rundir+'/Particles'
    path_to_input=rundir+"/input.dat"
    infile=open(path_to_input).readlines()
    for i in range(1,Nzones+1):
        for line in infile:
            if line.split("=")[0]==("computepart0%d:ngrains"%(i)):
                Nbins[i-1]=int(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:amin"%(i)):
                psizes_min[i-1]=float(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:amax"%(i)):
                psizes_max[i-1]=float(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:apow"%(i)):
                apows[i-1]=float(line.split("=")[1])

    for i in range(0,Nzones):
        psizes.append((psizes_min[i],psizes_max[i]))


    ############################################################
    # Creating lists of particle names
    fvC=str(fvC)
    fvSI=str(fvSI)
    case=[]
    for filename in os.listdir(folder):
        if fnmatch.fnmatch(filename,('*_f%s_f%s.fits.gz'%(fvSI,fvC))):
            case.append(filename)


    ############################################################
    # Define class 'Archivo'
    class Archivo:
        def __init__(self,filename):
            self.filename=filename

        def zone(self):
            return self.filename[8:12]

        def binnum(self):
            return self.filename[13:17]

        def f(self):
            hdulist=fits.open(folder+"/%s"%(self.filename))
            hdu=hdulist[0]
            hdr=hdu.header
            amin_bin=hdr["R_MIN"]
            amax_bin=hdr["R_MAX"]
            apow=hdr["R_POW"]
            z=int(self.zone())
            bins=Nbins[z-1]
            amin=psizes[z-1][0]
            amax=psizes[z-1][1]
            f_num=amax_bin**(-apow+4)-amin_bin**(-apow+4)
            f_den=amax**(-apow+4)-amin**(-apow+4)
            value=f_num/f_den
            return value

        def getEntry(self,i,j):
            hdulist=fits.open(folder+"/%s"%(self.filename))
            data=np.transpose(hdulist[0].data)
            value=data[i][j]
            return value

        def getWavelength(self):
            hdulist=fits.open(folder+"/%s"%(self.filename))
            data=np.transpose(hdulist[0].data)
            value=np.reshape(data[:,0:1],data.shape[0])
            return value

    ############################################################
    # Create wavelength array
    p=Archivo(case[0])
    wl=p.getWavelength()


    ############################################################
    # Create multidimensional matrices
    ext_matrix=np.zeros((Nzones,len(wl),int(max(Nbins))))
    abso_matrix=np.zeros((Nzones,len(wl),int(max(Nbins))))
    sca_matrix=np.zeros((Nzones,len(wl),int(max(Nbins))))


    ############################################################
    # Filling matrices
    Nbin_index=np.zeros(Nzones)
    for particle in case:
        p=Archivo(particle)
        zone_index=int(p.zone())
        j=int(Nbin_index[zone_index-1])
        for i in range(0,len(wl)):
            ext_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,1)
            abso_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,2)
            sca_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,3)
        Nbin_index[zone_index-1]+=1


    ############################################################
    # Adding entries
    sum_ext=np.zeros((Nzones,len(wl)))
    sum_abso=np.zeros((Nzones,len(wl)))
    sum_sca=np.zeros((Nzones,len(wl)))

    for j in range(0,Nzones):
        for i in range(0,len(wl)):
            sum_ext[j][i]=np.sum(ext_matrix[j][i])
            sum_abso[j][i]=np.sum(abso_matrix[j][i])
            sum_sca[j][i]=np.sum(sca_matrix[j][i])


    ############################################################
    # Creating HDU's 
    hdu_ext=np.zeros((Nzones+1,len(wl)))
    hdu_abso=np.zeros((Nzones+1,len(wl)))
    hdu_sca=np.zeros((Nzones+1,len(wl)))

    hdu_ext[0]=wl
    hdu_abso[0]=wl
    hdu_sca[0]=wl

    for j in range(0,Nzones):
        hdu_ext[j+1]=sum_ext[j]
        hdu_abso[j+1]=sum_abso[j]
        hdu_sca[j+1]=sum_sca[j]


    hdu_ext=fits.PrimaryHDU(np.transpose(hdu_ext))
    hdu_abso=fits.PrimaryHDU(np.transpose(hdu_abso))
    hdu_sca=fits.PrimaryHDU(np.transpose(hdu_sca))
    hdu_ext.writeto('ext.fits')
    hdu_abso.writeto('abs.fits')
    hdu_sca.writeto('sca.fits')

    return None



def plot_opacities(fits_ext,fits_abs,fits_sca):

    ############################################################
    #
    # Plot opacities extracted with "read_opacities()". 
    # Parameters: fits_ext, fits_abs, fits_sca [str]
    # containing the names of the files generated by 
    # "read_opacities()" 
    #
    ############################################################
    
    hdulist_ext=fits.open(fits_ext) 
    hdulist_abso=fits.open(fits_abs)
    hdulist_sca=fits.open(fits_sca)

    ext=hdulist_ext[0].data
    abso=hdulist_abso[0].data
    sca=hdulist_sca[0].data

    Nzones=ext.shape[1]-1

    ############################################################
    # Plotting
    fig=plt.figure(figsize=(12,5))
    gs=gridspec.GridSpec(1,Nzones,hspace=0.0)
    for i in range(0,Nzones):
        ax=plt.subplot(gs[0,i])
        ax.plot(ext[:,0:1],ext[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{ext}}$")
        ax.plot(abso[:,0:1],abso[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{abs}}$")
        ax.plot(sca[:,0:1],sca[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{sca}}$")
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=12,loc="lower left")
        ax.set_xlabel(r"$\lambda \, (\mu m)$")
        ax.set_ylim(5e-4,3e5)
        if i==0:
            ax.set_ylabel(r"Dust opacities (cm$^2$/g(dust))")
    plt.show()

    
