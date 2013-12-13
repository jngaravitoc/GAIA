import matplotlib.pyplot as plt
import numpy as np

dataO = np.loadtxt('MetalictiesSPH.dat') # Data original de Daisuke
ages= np.loadtxt('age.dat')

dataE = np.genfromtxt('ErMetalicitiesSPH.dat') # data con errores


Xreal = dataE[:, 17]
Yreal = dataE[:, 18]
Zreal = dataE[:, 19]

Xobs = dataE[:, 6]
Yobs = dataE[:, 7]
Zobs = dataE[:, 8]

#age = dataO[:, ]

G = dataE[:, 12]
GR = dataE[:, 13]

Rreal = np.sqrt(Xreal**2 + Yreal**2)
Robs = np.sqrt(Xobs**2 + Yobs**2)

ME = dataE[:, 14]
Mreal = dataE[:, 15]
Msd = dataE[:, 16]

res = np.where(Robs < 30.0)

ME_s = ME[res]
Mreal_s = Mreal[res]
Msd_s = Msd[res]
ages_s = ages[res]
 

Xobs_s = Xobs[res]
Yobs_s = Yobs[res]
Zobs_s = Zobs[res]

Xreal_s = Xreal[res]
Yreal_s = Yreal[res]
Zreal_s = Zreal[res]

Rreal_s = Rreal[res]
Robs_s = Robs[res]

plt.figure(1, figsize=(9.5, 9))
plt.hist2d(Rreal_s, Robs_s, bins=100)
p = np.polyfit(Rreal_s, Robs_s, 1)
plt.plot(Rreal, Rreal*p[0]+p[1], linewidth=2.5, color='r')
plt.xlim([0, 15])
plt.ylim([0, 20])
plt.xlabel('$\mathrm{R_{real}}$', fontsize=28)
plt.ylabel('$\mathrm{R_{obs}}$', fontsize=28)
plt.text(2, 13, '$\mathrm{Robs = Rreal\cdot 0.89 + 0.74}$', fontsize=28, color='r')
plt.savefig('rr.png')
plt.close()

print 'poly', p[0], p[1]

zselect = np.where(np.abs(Zobs_s <0.1))
zselectr = np.where(np.abs(Zreal_s <0.1))
RrealZ = Rreal_s[zselectr]
RobsZ = Robs_s[zselect]
MEZ = ME_s[zselect]
MrealZ = Mreal_s[zselectr]

z2select = np.where((np.abs(Zobs_s >0.1))&(Zobs_s<0.3))
z2selectr = np.where((np.abs(Zreal_s >0.1))&(Zreal_s<0.3))
RrealZ2 = Rreal_s[z2selectr]
RobsZ2 = Robs_s[z2select]
MEZ2 = ME_s[z2select]
MrealZ2 = Mreal_s[z2selectr]


plt.figure(1, figsize=(9.5, 9))
plt.scatter(RobsZ, MEZ)
plt.scatter(RrealZ, MrealZ, c='r')
plt.savefig('metalicitiZ.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.scatter(RobsZ2, MEZ2)
plt.scatter(RrealZ2, MrealZ2, c='r')
plt.savefig('metalicitiZ2.png')
plt.close()


print np.amax(Robs), np.amax(Xobs), np.amax(Yobs)
print np.amax(Rreal)
print np.average(Robs)

nbin=15
rmax=15.0

plt.figure(1, figsize=(9.5, 9))
plt.hexbin(Xreal_s, Yreal_s, Mreal_s)
plt.ylim([-15, 15])
plt.xlim([-15, 20])
plt.xlabel('$\mathrm{X(Kpc)}$', fontsize=28)
plt.ylabel('$\mathrm{Y(Kpc)}$', fontsize=28)
plt.colorbar()
plt.savefig('realMet.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.hexbin(Xobs_s, Yobs_s, ME_s, vmin=-0.8 , vmax=0.2)
plt.ylim([-15, 15])
plt.xlim([-15, 20])
plt.xlabel('$\mathrm{X(Kpc)}$', fontsize=28)
plt.ylabel('$\mathrm{Y(Kpc)}$', fontsize=28)
plt.colorbar()
plt.savefig('obsMet.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.hexbin(Xobs_s, Yobs_s, ages_s)
plt.colorbar()
plt.ylim([-12, 12])
plt.xlim([-15, 20])
plt.savefig('obsAge.png')
plt.close()

np_r=np.histogram(Robs_s,nbin,(0.0,15.0))[0]
# mmean metallicity
MmE=np.histogram(Robs_s,nbin,(0.0,15.0),weights=ME_s)[0]/np_r
# square mean
Mm2E=np.histogram(Robs_s,nbin,(0.0,15.0),weights=(ME_s**2))[0]/np_r
#Metalicity disperison 
MsigE = (Mm2E - MmE**2)**0.5


np_e=np.histogram(Rreal_s,nbin,(0.0,15.0))[0]
# mmean metallicity
Mm=np.histogram(Rreal_s,nbin,(0.0,rmax),weights=Mreal_s)[0]/np_e
# square mean
Mm2=np.histogram(Rreal_s,nbin,(0.0,rmax),weights=(Mreal_s**2))[0]/np_e
#Metalicity disperison 
Msig = (Mm2 - Mm**2)**0.5



# Haciendo los anillos para los datos sin errores
r = Rreal_s 
r1 = np.where(r<1.0)
r2 = np.where((r>1.0) & (r<2.0))
r3 = np.where((r>2.0) & (r<3.0))
r4 = np.where((r>3.0) & (r<4.0))
r5 = np.where((r>4.0) & (r<5.0))
r6 = np.where((r>5.0) & (r<6.0))
r7 = np.where((r>6.0) & (r<7.0))
r8 = np.where((r>7.0) & (r<8.0))
r9 = np.where((r>8.0) & (r<9.0))
r10 = np.where((r>9.0) & (r<10.0))
r11= np.where((r>10.0) & (r<11.0))
r12= np.where((r>11.0) & (r<12.0))
r13= np.where((r>12.0) & (r<13.0))
r14= np.where((r>13.0) & (r<14.0))
r15= np.where((r>14.0) & (r<15.0))

N_r = [len(ME_s[r1]), len(ME_s[r2]), len(ME_s[r3]), len(ME_s[r4]), len(ME_s[r5]), len(ME_s[r6]), len(ME_s[r7]), len(ME_s[r8]), len(ME_s[r9]), len(ME_s[r10]), len(ME_s[r11]), len(ME_s[r12]), len(ME_s[r13]), len(ME_s[r14]), len(ME_s[r15])]


rE = Robs_s

r1E = np.where(rE<1.0)
r2E = np.where((rE>1.0) & (rE<2.0))
r3E = np.where((rE>2.0) & (rE<3.0))
r4E = np.where((rE>3.0) & (rE<4.0))
r5E = np.where((rE>4.0) & (rE<5.0))
r6E = np.where((rE>5.0) & (rE<6.0))
r7E = np.where((rE>6.0) & (rE<7.0))
r8E = np.where((rE>7.0) & (rE<8.0))
r9E = np.where((rE>8.0) & (rE<9.0))
r10E = np.where((rE>9.0) & (rE<10.0))
r11E= np.where((rE>10.0) & (rE<11.0))
r12E= np.where((rE>11.0) & (rE<12.0))
r13E= np.where((rE>12.0) & (rE<13.0))
r14E= np.where((rE>13.0) & (rE<14.0))
r15E= np.where((rE>14.0) & (rE<15.0))

N_E = [len(ME_s[r1E]), len(ME_s[r2E]), len(ME_s[r3E]), len(ME_s[r4E]), len(ME_s[r5E]), len(ME_s[r6E]), len(ME_s[r7E]), len(ME_s[r8E]), len(ME_s[r9E]), len(ME_s[r10E]), len(ME_s[r11E]), len(ME_s[r12E]), len(ME_s[r13E]), len(ME_s[r14E]), len(ME_s[r15E])]


# Magnitudes promedio en cada anillo sin errores
M0=Mreal_s

M1o = np.average(M0[r1])
M2o = np.average(M0[r2])
M3o = np.average(M0[r3])
M4o = np.average(M0[r4])
M5o = np.average(M0[r5])
M6o = np.average(M0[r6])
M7o = np.average(M0[r7])
M8o = np.average(M0[r8])
M9o = np.average(M0[r9])
M10o = np.average(M0[r10])
M11o = np.average(M0[r11])
M12o = np.average(M0[r12])
M13o = np.average(M0[r13])
M14o = np.average(M0[r14])
M15o = np.average(M0[r15])

ME = ME_s

M1E = np.average(ME[r1])
M2E = np.average(ME[r2])
M3E = np.average(ME[r3])
M4E = np.average(ME[r4])
M5E = np.average(ME[r5])
M6E = np.average(ME[r6])
M7E = np.average(ME[r7])
M8E = np.average(ME[r8])
M9E = np.average(ME[r9])
M10E = np.average(ME[r10])
M11E = np.average(ME[r11])
M12E = np.average(ME[r12])
M13E = np.average(ME[r13])
M14E = np.average(ME[r14])
M15E = np.average(ME[r15])

plt.hist(M0[r1], bins=40, normed=True, edgecolor='r', histtype='stepfilled', fill=False)
plt.hist(ME[r1], bins=40, normed=True, edgecolor='b', histtype='stepfilled', fill=False)
plt.title('$\mathrm{R_{gal}=[0-1Kpc]}$', fontsize=28)
plt.xlabel('$\mathrm{[Fe/H]}$', fontsize=28)
plt.savefig('hist1.png')
plt.close()


MTE = [M1E, M2E, M3E, M4E, M5E, M6E, M7E, M8E, M9E, M10E, M11E, M12E, M13E, M14E, M15E]
MTO = [M1o, M2o, M3o, M4o, M5o, M6o, M7o, M8o, M9o, M10o, M11o, M12o, M13o, M14o, M15o]

E = [np.abs(M1o-M1E), np.abs(M2o-M2E), np.abs(M3o-M3E), np.abs(M4o-M4E), np.abs(M5o-M5E),np.abs(M6o-M6E),np.abs(M7o-M7E),np.abs(M8o-M8E),np.abs(M9o-M9E), np.abs(M10o-M10E),np.abs(M11o-M11E),np.abs(M12o-M12E),np.abs(M13o-M13E) ,np.abs(M14o-M14E),np.abs(M15o-M15E)]

rO = ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

t = np.linspace(0, 15, 15)
plt.figure(1, figsize=(9.5, 9))
plt.plot(rO, MTO, c='r', marker='^')
plt.plot(rO, MTE, marker='o', c='b')
plt.errorbar(rO, MTE, yerr = MsigE, fmt='ob')
plt.errorbar(rO, MTO, yerr = Msig, fmt='or')
p = np.polyfit(rO, MTO, 1)
#plt.plot(t,p[0]*t+p[1])
plt.xlabel(r'$\mathrm{R_{gal}}$', fontsize=28)
plt.ylabel(r'$\mathrm{[Fe/H]}$', fontsize=28)
plt.legend(['No errors','With errors'])
#plt.show()
plt.savefig('toericavsexperimental.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.plot(rO, N_r, c ='k', linewidth=1.5)
#plt.plot(rO, N_E, c='red', linewidth=1.5)
plt.xlabel('$\mathrm{R_{gal}}$', fontsize=28)
plt.ylabel('$\mathrm{N_{stars}}$', fontsize = 28)
#plt.legend(['real', 'errors'])
plt.savefig('histRealandobs.png')
plt.close()

Gt =np.linspace(0, 18, 36)
ape= 0.81192 - 0.12226 *Gt + 0.0056669 *Gt**2
apem = 0.71709 - 0.115 * Gt + 0.0050649*Gt**2

plt.figure(1, figsize=(9.5, 9))
plt.plot(Gt, apem, c='r')
plt.plot(Gt, ape, c = 'g')
plt.legend(['$\mathrm{Internal}$ $\mathrm{Errors}$', '$\mathrm{External}$ $\mathrm{Errors}$'])
plt.scatter(G, Msd)
plt.ylabel('$\mathrm{\sigma_{Fe/H}}$', fontsize=28)
plt.xlabel('$\mathrm{G}$', fontsize=28)
plt.savefig('GvsDelta.png')
plt.close()

plt.plot(rO, N_E)
plt.savefig('histObs.png')
plt.close()
"""
plt.plot(rO, E, marker='o')
plt.xlabel(r'$\mathrm{R_{gal}}$', fontsize=28)
plt.ylabel(r'$\mathrm{Log(Z)}$', fontsize=28)
plt.legend(['Errors'])
#plt.show()
plt.savefig('errors.png')
plt.close()
"""
print len(Rreal)
print len(Rreal_s)

#plt.scatter(Rreal_s, Mreal_s)
#lt.savefig('experimental.png')
#plt.xlim([0,15])
#plt.ylim([0,15])
#plt.savefig('RealZ.png')

#plt.hist2d(Robs_s, ME_s, bins=160)
x_teo = np.linspace(0, 15, 15)
y_teo = -0.05*x_teo + 0.2

plt.figure(1, figsize=([9.5, 9]))
plt.scatter(Robs_s, ME_s,s=5,  c='b', alpha=0.5)
plt.scatter(Rreal_s, Mreal_s,s=5, c = 'r', alpha=0.5)
plt.plot(x_teo, y_teo,  linewidth=3.5, c='g')
plt.plot(t,p[0]*t+p[1], c='k', linewidth=3.5)
plt.legend(['Inital condition', 'With errors', 'No errors', 'Observed Regresion'])
#lt.axhline(0.0, 15, linewidth=3.5, c='m')
#plt.plot(t,p[0]*t+p[1])
plt.xlim([0, 15])
plt.ylim([-4.0, 4.0])
plt.xlabel(r'$\mathrm{R_{gal}}$', fontsize=28)
plt.ylabel(r'$\mathrm{[Fe/H]}$', fontsize=28)
#plt.legend(['With errors', 'No errors', 'Initial Condition'])
plt.scatter(8, 0, marker='*', s=190, c='y')
plt.grid()
plt.savefig('ObsZandreal.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.scatter((np.sqrt((Xobs_s + 8.5)**2 + Yobs_s**2)), np.abs(ME_s-Mreal_s), s=5, c='b', alpha=0.5)
plt.xlabel('$\mathrm{D_{\odot}}$', fontsize=28)
plt.ylabel('$\mathrm{|(Fe/H)_{o} -((Fe/H)_{r})|}$', fontsize = 28)
plt.xlim([0, 40])
plt.ylim([0, 1.4])
plt.savefig('ErrosvsDsol.png')
plt.close()

plt.figure(1, figsize=(9.5, 9))
plt.hist2d((np.sqrt((Xobs_s + 8.5)**2 + Yobs_s**2)), np.abs(ME_s-Mreal_s), bins=90)
#plt.xlabel('$\mathrm{D_{\odot}}$', fontsize=28)
#plt.ylabel('$\mathrm{|(Fe/H)_{o} -((Fe/H)_{r})|}$', fontsize = 28)
plt.xlim([0, 20])
plt.ylim([0, 1.2])
plt.savefig('HistErrosvsDsol.png')
plt.close()


plt.figure(1,figsize=(9.5, 9))
#plt.scatter(np.sqrt((Xobs_s)**2 + Yobs_s**2), np.abs(Robs_s - Rreal_s), c='r', alpha=0.1)
plt.scatter(np.sqrt((Xobs_s + 8.5)**2 + Yobs_s**2), np.abs(Robs_s - Rreal_s),s=0.2, alpha=0.5)
#plt.scatter(np.sqrt((Xobs_s)**2 + Yobs_s**2), np.abs(Robs_s - Rreal_s), c='r')
plt.xlabel('$\mathrm{D_{\odot}}$', fontsize=28)
plt.ylabel('$\mathrm{|R_{o} - R_{r}|}$', fontsize = 28)
plt.xlim([0, 40])
plt.ylim([0, 30])
plt.savefig('radiodifvsDsol.png')
plt.close()


# Grafica 

"""
plt.hist2d(xO, yO,  bins=80)
plt.show()
"""
