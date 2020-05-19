import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch


############################################################
# Inputs
Nzones=3
run_name="run_166"
fvSI=0.73
fvC=0.27


############################################################
# Initializing stuff
Nbins=np.zeros(Nzones)
apows=np.zeros(Nzones)
psizes_min=np.zeros(Nzones)
psizes_max=np.zeros(Nzones)
psizes=[]


############################################################
# Working dir
directory=("/data/users/bportilla/runs/final_runs/%s/"%(run_name))
folder=directory+'Particles'
path_to_input=directory+"input.dat"
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
# Creating lists for both cases
fvC=str(fvC)
fvSI=str(fvSI)
case2=[]
for filename in os.listdir(folder):
    if fnmatch.fnmatch(filename,('*_f%s_f%s.fits.gz'%(fvSI,fvC))):
        case2.append(filename)


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
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/%s/Particles/%s"%(run_name,self.filename))
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
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/%s/Particles/%s"%(run_name,self.filename))
        data=np.transpose(hdulist[0].data)
        value=data[i][j]
        return value
    
    def getWavelength(self):
        hdulist=fits.open("/data/users/bportilla/runs/final_runs/%s/Particles/%s"%(run_name,self.filename))
        data=np.transpose(hdulist[0].data)
        value=np.reshape(data[:,0:1],data.shape[0])
        return value


############################################################
# Create wavelength array
p2=Archivo(case2[0])
wl2=p2.getWavelength()



############################################################
# Create matrices to store values of Ai,Bi and Ci for case 1
# and case 2

ext_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))
abso_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))
sca_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))



############################################################
# Filling matrices
Nbin_index=np.zeros(Nzones)
for particle in case2:
    p=Archivo(particle)
    zone_index=int(p.zone())
    j=int(Nbin_index[zone_index-1])
    for i in range(0,len(wl2)):
        ext_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,1)
        abso_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,2)
        sca_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,3)
    Nbin_index[zone_index-1]+=1
    
  

############################################################
# Adding entries
sum_ext=np.zeros((Nzones,len(wl2)))
sum_abso=np.zeros((Nzones,len(wl2)))
sum_sca=np.zeros((Nzones,len(wl2)))

for j in range(0,Nzones):
    for i in range(0,len(wl2)):
        sum_ext[j][i]=np.sum(ext_matrix[j][i])
        sum_abso[j][i]=np.sum(abso_matrix[j][i])
        sum_sca[j][i]=np.sum(sca_matrix[j][i])


############################################################
# Creating HDU's for case 2
hdu_ext=np.zeros((Nzones+1,len(wl2)))
hdu_abso=np.zeros((Nzones+1,len(wl2)))
hdu_sca=np.zeros((Nzones+1,len(wl2)))

hdu_ext[0]=wl2
hdu_abso[0]=wl2
hdu_sca[0]=wl2

for j in range(0,Nzones):
    hdu_ext[j+1]=sum_ext[j]
    hdu_abso[j+1]=sum_abso[j]
    hdu_sca[j+1]=sum_sca[j]

print(np.transpose(hdu_ext).shape)
print(np.transpose(hdu_ext))

hdu_ext=fits.PrimaryHDU(np.transpose(hdu_ext))
hdu_abso=fits.PrimaryHDU(np.transpose(hdu_abso))
hdu_sca=fits.PrimaryHDU(np.transpose(hdu_sca))
hdu_ext.writeto('ext.fits')
hdu_abso.writeto('abs.fits')
hdu_sca.writeto('sca.fits')




