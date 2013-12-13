
# to analyse rotation and radial mean velocity and dispersion
# from error added data from UB code

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants as const
from FortranFile import FortranFile
import struct

# input parameters
# for radial velocity profile
# number of bin
nbin=12
# number of rmax (kpc)
rmax=12.0
# magnitude limit
vmaglim=16.5
# assumed rotation velocity at 8 kpc
vrotsun=220.0
# position of the sun
xsun=-8.0

# input file name
inputfilei='../ubgaiaerr/gaiaei-out.dat'
inputfiled='../ubgaiaerr/gaiaed-out.dat'

# reading binary data
f=FortranFile(inputfilei)
nset=f.readInts()
print ' input Number of stars =',nset
f.close()
f=FortranFile(inputfiled)
i=0
rdata=np.zeros((nset,25))
while i < nset:
  rdata[i,:]=f.readReals('d')
  i+=1
f.close()


# parallax (mas) -> arcsec
paradet=rdata[:,8]/1000.0
# only chose positive parallax and < 16 kpc away
sindx=np.where(paradet > 1.0/(16.0*1.0e3))

# counting using parallax
parado=paradet[sindx]
nse=len(parado)
print ' Number of stars (D<16 kpc)=',nse

# set the data  ...o (error added), ...e (error)
# alpha, delta (radian)
alpha=np.reshape(rdata[sindx,0],nse)
delta=np.reshape(rdata[sindx,1],nse)
alphao=np.reshape(rdata[sindx,6],nse)
deltao=np.reshape(rdata[sindx,7],nse)
alphae=np.reshape(rdata[sindx,12],nse)
deltae=np.reshape(rdata[sindx,13],nse)
# parallax mas -> (arcsec)
parad=np.reshape(rdata[sindx,2],nse)*0.001
parade=np.reshape(rdata[sindx,14],nse)*0.001
# proper motion mas/yr -> (arcsec/yr)
mua=np.reshape(rdata[sindx,3],nse)*0.001
mud=np.reshape(rdata[sindx,4],nse)*0.001
muao=np.reshape(rdata[sindx,9],nse)*0.001
mudo=np.reshape(rdata[sindx,10],nse)*0.001
muae=np.reshape(rdata[sindx,15],nse)*0.001
mude=np.reshape(rdata[sindx,16],nse)*0.001
# radial velocity (km/s)
vrad=np.reshape(rdata[sindx,5],nse)
vrado=np.reshape(rdata[sindx,11],nse)
vrade=np.reshape(rdata[sindx,17],nse)
# Gmag
gmag=np.reshape(rdata[sindx,18],nse)
# G(BP-RP)
gcol=np.reshape(rdata[sindx,19],nse)
# [Fe/H]
feh=np.reshape(rdata[sindx,20],nse)
# Av
ave=np.reshape(rdata[sindx,21],nse)
#
vmag=np.reshape(rdata[sindx,22],nse)
# true dist
dist=np.reshape(rdata[sindx,23],nse)
# GRVS
grvs=np.reshape(rdata[sindx,24],nse)


# below ignoring the errror, but just taking the face values

# get Cartician position in equatorial coordidate
# distance in kpc
disto=1.0/parado/1000.0
# x,y,z position from the sun's position
rxy=disto*np.cos(deltao)
px=rxy*np.cos(alphao)
py=rxy*np.sin(alphao)
pz=disto*np.sin(deltao)
# vx,vy,vz
# arcsec/yr -> km/s
valpha=4.74*muao/parado
vdelta=4.74*mudo/parado
vxy=vrado*np.cos(deltao)-vdelta*np.sin(deltao)
vx=vxy*np.cos(alphao)-valpha*np.sin(alphao)
vy=vxy*np.sin(alphao)+valpha*np.cos(alphao)
vz=vrado*np.sin(deltao)+vdelta*np.cos(deltao)

#f=open('testgaein-part.dat','w')
#i=0
#while i < nse:
# print >>f, "%f %f %f %f %f %f %f %f %f" %(px[i],py[i],pz[i]
#  ,vx[i],vy[i],vz[i],alphao[i],deltao[i],disto[i])
# i+=1
#f.close()


# convert to galactic coordinate
# 3 axes in Galactic coordinate
# North Celestial Pole (z-axis)
lzeq=np.radians(122.93193212)
bzeq=np.radians(27.12835496)
# RA,DEC=0,0 (x-axis)
lxeq=np.radians(96.33723825)
bxeq=np.radians(-60.18853909)
#  RA,DEC=90,0 (y-axis)
lyeq=np.radians(206.98916373)
byeq=np.radians(-11.42442440)
# transformation matrix (though it is array in python)
tmateqga=np.array([
  [np.cos(lxeq)*np.cos(bxeq),np.sin(lxeq)*np.cos(bxeq),np.sin(bxeq)],
  [np.cos(lyeq)*np.cos(byeq),np.sin(lyeq)*np.cos(byeq),np.sin(byeq)],
  [np.cos(lzeq)*np.cos(bzeq),np.sin(lzeq)*np.cos(bzeq),np.sin(bzeq)]
])
print 'tmat gal->eq=',tmateqga
# take inverse
tmatgaeq=np.linalg.inv(tmateqga)
print 'tmat eq->gal=',tmatgaeq

print 'tmat gal->eq= %12.9f %12.9f %12.9f' % (tmateqga[0,0],tmateqga[0,1],tmateqga[0,2])
print '              %12.9f %12.9f %12.9f' % (tmateqga[1,0],tmateqga[1,1],tmateqga[1,2])
print '              %12.9f %12.9f %12.9f' % (tmateqga[2,0],tmateqga[2,1],tmateqga[2,2])

print 'tmat eq->gal= %12.9f %12.9f %12.9f' % (tmatgaeq[0,0],tmatgaeq[0,1],tmatgaeq[0,2])
print '              %12.9f %12.9f %12.9f' % (tmatgaeq[1,0],tmatgaeq[1,1],tmatgaeq[1,2])
print '              %12.9f %12.9f %12.9f' % (tmatgaeq[2,0],tmatgaeq[2,1],tmatgaeq[2,2])



# transfer the coordinate: x-y plane: Cerestial equator -> Galactic plane
# position
# position in Galactic cartesian coordinate
poseqs=np.vstack([px,py,pz])
# transfer to equatorial cartesian coordinate
posgas=np.dot(tmatgaeq,poseqs)

# velocity
# velocity in Galactic cartesian coordinate
veleqs=np.vstack([vx,vy,vz])
# transfer to equatorial cartesian coordinate
velgas=np.dot(tmatgaeq,veleqs)

# reset position and velocity
px=posgas[0,:]+xsun
py=posgas[1,:]
pz=posgas[2,:]
vx=velgas[0,:]
vy=velgas[1,:]+vrotsun
vz=velgas[2,:]

rad=np.sqrt(np.power(px,2)+np.power(py,2))

print ' rmax=',np.max(rad)

# rotation velocity 
vrot=(vx*py-vy*px)/rad
# radial velocity 
vrad=(vx*px+vy*py)/rad

# output ascii data for test
#f=open('testgae-part.dat','w')
#i=0
#while i < nse:
# print >>f, "%f %f %f %f %f %f %f %f %f" %(px[i],py[i],pz[i]
#  ,vx[i],vy[i],vz[i],vrot[i],vrad[i],rad[i])
# i+=1
#f.close()

# velocity mean and dispersion
np_r=np.histogram(rad,nbin,(0.0,rmax))[0]
# mean radius in each bin
rmean_r=np.histogram(rad,nbin,(0.0,rmax),weights=rad)[0]/np_r

# Rotation velocity
# mean
vrotm_r=np.histogram(rad,nbin,(0.0,rmax),weights=vrot)[0]/np_r
# square mean
vrotm2_r=np.histogram(rad,nbin,(0.0,rmax),weights=(vrot**2))[0]/np_r
# velocity dispersion 
vrotsig_r=(vrotm2_r-vrotm_r**2)**0.5

# Radial velocity
# mean
vradm_r=np.histogram(rad,nbin,(0.0,rmax),weights=vrad)[0]/np_r
# square mean
vradm2_r=np.histogram(rad,nbin,(0.0,rmax),weights=(vrad**2))[0]/np_r
# velocity dispersion
vradsig_r=(vradm2_r-vradm_r**2)**0.5

# output ascii data
f=open('velerr_r.dat','w')
i=0
while i < nbin:
 print >>f, "%f %f %f %f %f %d" %(rmean_r[i],vrotm_r[i],vrotsig_r[i]
  ,vradm_r[i],vradsig_r[i],np_r[i])
 i+=1
f.close()

plt.subplot(111)
# labes
plt.xlabel("R (kpc)",fontsize=18,fontname="serif")
plt.ylabel(r"$\rm V_{rot}$",fontsize=18,fontname="serif",style="normal")
# hexbin plot
plt.hexbin(rad,vrot,bins='log',gridsize=300,cmap=cm.jet)
# plot mean
plt.errorbar(rmean_r,vrotm_r,yerr=vrotsig_r,fmt='ow')
#plt.axis([pxmin,pxmax,pymin,pymax])
plt.axis([0.0,rmax,0.0,300.0])
#cb=plt.colorbar()
plt.show()
plt.savefig('velerr_r.png')
