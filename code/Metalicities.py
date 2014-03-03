
# replace N-body star particle to RC star (Mv=1.2714, M_G=1, .V-I=1.0)
# generate UB Gaia error fortran code input file.
# for stars their V (including extinction) <16.5
# 

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from scipy import constants as const
from FortranFile import FortranFile
import struct

# output option 0: off, 1: on
# output UB Gaia error code input binary data
flagub=1
# plot 
flagplot=1

# input filename
# ~1/100 data
#inputfile='../simdata/disk1s.dat'
# number of particle
nset=42075
# full data
inputfile='../data/disk1.dat'
# number of particle
nset=4207480

# magnitude limit
vmaglim=16.5
# assumed rotation velocity at 8 kpc
vrotsun=220.0
# assumed M_v of RC, V-I colour
mvrc=1.2714
virc=1.0
# solar metallicity
zsun=np.float(0.019)

f=FortranFile(inputfile)
i=0
rdata=np.zeros((nset,10))
while i < nset:
  rdata[i,:]=f.readReals('d')
  i+=1
f.close()

print ' Number of stars=',nset

# set date
# position (kpc)
xs=rdata[:,0]
ys=rdata[:,1]
zs=rdata[:,2]
# velocity (km/s)
vxs=rdata[:,3]
vys=rdata[:,4]
vzs=rdata[:,5]
# particle mass (Msun)
ms=rdata[:,6]
# age (Gyr)
ages=rdata[:,7]
# metallicity  (Z/Zsun)
mets=rdata[:,8]/zsun
# Av
avs=rdata[:,9]*5.0

# distance in kpc
rads=np.sqrt(np.power(xs,2)+np.power(ys,2)+np.power(zs,2))

# apparent magnitude
vmag=mvrc-5.0*(1.0-np.log10(rads*1000.0))
imag=vmag-virc
# I band extinction (from Galaxia, assuminv R_V=3.1)
ais=0.594*avs
# add extinction
vmage=vmag+avs
image=imag+ais
vicole=vmage-image

# selecting apparent magnitude < vmaglim
sindx=np.where(vmage < vmaglim)
ns=np.size(sindx)
print 'Ns(V<',vmaglim,')=',ns
xos=xs[sindx]
yos=ys[sindx]
zos=zs[sindx]
vxos=vxs[sindx]
vyos=vys[sindx]
vzos=vzs[sindx]

# ages
ageos=ages[sindx]
# [Fe/H] = logZ/Zsun
fehos=np.log10(mets[sindx])

#ages2 = ages[sindx] is already defined
#mets2 = mets[sindx]
avos=avs[sindx]
vmageos=vmage[sindx]
vicoleos=vicole[sindx]
"""
zselect = np.where((np.abs(zos)>0.0) & (np.abs(zos)<0.3))

nsz = np.size(zselect)

xosz=xs[zselect]
yosz=ys[zselect]
zosz=zs[zselect]
vxosz=vxs[zselect]
vyosz=vys[zselect]
vzosz=vzs[zselect]
ageosz=ages[zselect]
avosz=avs[zselect]
fehosz=np.log10(mets[zselect])
# Galactic coordinate
# distance in the plane in kpc
radxyos=np.sqrt(np.power(xos,2)+np.power(yos,2))
# distance in 3D in kpc
"""
"""
rados=np.sqrt(np.power(xos,2)+np.power(yos,2)+np.power(zos,2))
# longitude
los=np.degrees(np.arccos(xos/radxyos))
i=0
while i<ns:
 if yos[i] < 0.0:
   los[i]=360.0-los[i]
 i+=1
# latitude
bos=np.degrees(np.arcsin(zos/rados))
galos=np.vstack([los,bos])

# transfer the coordinate: x-y plane: Galactic plane -> Cerestial equator

# setting up a transfer matrix
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

# position
# position in Galactic cartesian coordinate
posgaos=np.vstack([xos,yos,zos])
# transfer to equatorial cartesian coordinate
poseqos=np.dot(tmateqga,posgaos)

# velocity
# velocity in Galactic cartesian coordinate
velgaos=np.vstack([vxos,vyos,vzos])
# transfer to equatorial cartesian coordinate
veleqos=np.dot(tmateqga,velgaos)

# R.A. Dec. 
radxyos=np.sqrt(np.power(poseqos[0,:],2)+np.power(poseqos[1,:],2))
rados=np.sqrt(np.power(poseqos[0,:],2)+np.power(poseqos[1,:],2)
  +np.power(poseqos[2,:],2))
# R.A.
alpos=np.degrees(np.arccos(poseqos[0,:]/radxyos))
i=0
while i<ns:
 if poseqos[1,i] < 0.0:
   alpos[i]=360.0-alpos[i]
 i+=1
delos=np.degrees(np.arcsin(poseqos[2,:]/rados))
equos=np.vstack([alpos,delos])

# radial velocity
vrados=(veleqos[0,:]*poseqos[0,:]+veleqos[1,:]*poseqos[1,:]
  +veleqos[2,:]*poseqos[2,:])/rados
# proper motion
# mu_alpha (km/s)
valpos=(-veleqos[0,:]*poseqos[1,:]+veleqos[1,:]*poseqos[0,:])/radxyos
# mu_delta (km/s)
# vrad in x-y plane
vradxyos=(veleqos[0,:]*poseqos[0,:]+veleqos[1,:]*poseqos[1,:])/radxyos
vdelos=(-vradxyos*poseqos[2,:]+veleqos[2,:]*radxyos)/rados

# changing the unit from km/s to mas/yr
# km/s -> arcsec/y   vt=4.74 mu d(pc)
valpos=1000.0*valpos/4.74/(rados*1000.0)
vdelos=1000.0*vdelos/4.74/(rados*1000.0)
"""
# test ASCII outputt
f=open('../data/MetalictiesSPH.dat','w')
i=0
print >>f, ("# X  Y  Z  VX  VY VZ AVOS FEHOS AGES")
while i < ns:
 print >>f, "%f %f %f %f %f %f %f %f %f" %(xos[i] - 8.0 , yos[i],zos[i],vxos[i] ,vyos[i],vzos[i], avos[i], fehos[i], ages[i])
 i+=1
f.close()

"""
plt.figure(1, figsize=(9.5, 9))
plt.hist2d(xos - 8.0, yos, bins=290)
plt.title('$\mathrm{MW}$ $\mathrm{RC}$ $\mathrm{Densities}$', fontsize=28)
plt.ylabel('$\mathrm{Y(Kpc)}}$', fontsize=28)
plt.xlabel('$\mathrm{X(Kpc)}}$', fontsize=28)
plt.scatter(-8.5, 0, s=190,c= 'y', marker= '*')
#plt.xlim([-14, 14])
#plt.ylim([-14, 14])
plt.savefig('GalaxyRC.png')
plt.close()

#plt.plot(rO, N_r)
#plt.plot(rO, N_E, c='r')
plt.figure(1, figsize=(9.5, 9))
plt.hexbin(xos-8.0, yos,fehos )
plt.colorbar()
plt.ylabel('$\mathrm{Y(Kpc)}}$', fontsize=28)
plt.xlabel('$\mathrm{X(Kpc)}}$', fontsize=28)
plt.title('$\mathrm{Metallicity}$ $\mathrm{distribution}$', fontsize=28)
plt.savefig('metsdistribution.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.hexbin(xos-8.0, yos, ageos)
plt.colorbar()
plt.ylabel('$\mathrm{Y(Kpc)}}$', fontsize=28)
plt.xlabel('$\mathrm{X(Kpc)}}$', fontsize=28)
plt.title('$\mathrm{Age}$ $\mathrm{distribution}$', fontsize=28)
plt.savefig('agedistribution.png')
plt.close()
"""

# distance in parsec
#diss=rads*1000.0
# binary data for UB Fortran code
if flagub:
  f=FortranFile('ubgaiaein.bin',mode='w')
  nsarr=np.reshape(ns,1)
  f.writeInts(nsarr)
  i=0
  while i < ns:
    staro=np.array([alpos[i],delos[i],rados[i],valpos[i],vdelos[i],vrados[i]
     ,fehos[i],avos[i],vicoleos[i],vmageos[i]], ages[i])
    f.writeReals(staro,prec='d')
    i+=1
  f.close()
"""
if flagplot:
# test plot
# top panels
  gs1=gridspec.GridSpec(1,3)
  gs1.update(left=0.05,right=0.975,bottom=0.5,top=0.95,hspace=0.5)
# x-y plot
  plt.subplot(gs1[0],aspect='equal')
  plt.hexbin(xos,yos,bins='log',gridsize=100,cmap=cm.jet)
  plt.axis([-5.0,15.0,-10.0,10.0])
  plt.xlabel("x-y",fontsize=12,fontname="serif")
  plt.ylabel("y",fontsize=12,fontname="serif")
# x-z plot
  plt.subplot(gs1[1],aspect='equal')
  plt.hexbin(xos,zos,bins='log',gridsize=100,cmap=cm.jet)
  plt.axis([-5.0,15.0,-10.0,10.0])
  plt.xlabel("x-z",fontsize=12,fontname="serif")
  plt.ylabel("z",fontsize=12,fontname="serif")
# Mv vs. V-I plotg
  plt.subplot(gs1[2],aspect=0.2)
# hexbin plot
  plt.hexbin(vicoleos,vmageos,bins='log',gridsize=200,cmap=cm.jet)
  plt.axis([1.0,2.5,16.5,6.0])
# labes
  plt.xlabel(r'$\rm V vs. V-I$',fontsize=12,fontname="serif")
  plt.ylabel(r'$\rm M_V$',fontsize=12,fontname="serif")

# bottom panel
  gs2=gridspec.GridSpec(1,2)
  gs2.update(left=0.1,right=0.975,bottom=0.05,top=0.4)
# Galactic coordinate
  plt.subplot(gs2[0],aspect='equal')
  plt.hexbin(los,bos,bins='log',gridsize=100,cmap=cm.jet)
  plt.axis([0.0,360.0,-90.0,90.0])
  plt.xlabel("l",fontsize=12,fontname="serif")
  plt.ylabel("b",fontsize=12,fontname="serif")
# equatorial coordinate
  plt.subplot(gs2[1],aspect='equal')
  #plt.scatter(np.sqrt(xos**2 + yos**2), fehos)
  plt.hexbin(equos[0,:],equos[1,:], bins='log', gridsize=100, cmap=cm.jet)
  #plt.axis([0.0,360.0,-90.0,90.0])
  plt.xlabel(r"$\alpha$",fontsize=12,fontname="serif")
  plt.ylabel(r"$\delta$",fontsize=12,fontname="serif")
#cb=plt.colorbar()
  #plt.show()
  plt.savefig('RCsamp.png')
"""
