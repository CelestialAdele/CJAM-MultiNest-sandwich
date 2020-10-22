import numpy as np
import matplotlib.pyplot as plt



x,y = np.loadtxt('xy.dat', dtype = 'float', usecols=(0,1), delimiter = ",", unpack = True)

vx,vy,vz,v2xx,v2yy,v2zz,v2xy,v2yz,v2xz,sigx,sigy,sigz = np.loadtxt('nc_new_moments.dat', dtype='float', usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack = True)

#UNITS: vx (mas/yr), vy (mas/yr), vz (km/s), v2xx (mas^2/yr^2), v2yy (mas^2/yr^2), v2zz(km^2/s^2), v2xy (mas^2,yr^2), v2xz (km*mas/yr*s), v2yz (km*mas/yr*s), x (arcsec), y (arcsec)


#convert x,y to arcmins

x = x*0.0167
y = y*0.0167

fig, ax2 = plt.subplots()

ax2.plot(x,y, 'ko', ms=3, alpha = 0.1)


#convert dispersions to (km/s)

v2xx = np.sqrt(v2xx)*4.74*4.6
v2yy = np.sqrt(v2yy)*4.74*4.6
v2zz = np.sqrt(v2zz)
v2xz = v2xz*4.74*4.6
v2xy = v2xy*4.74*4.6*4.74*4.6
v2xy = np.array(v2xy)

v2flags = np.zeros(len(v2xy))
	
for i in range(len(v2xy)):
	if v2xy[i] < 0.0:
		v2xy[i] = v2xy[i]*(-1.0)
		v2flags[i] = 1

v2xy = np.sqrt(v2xy)

for i in range(len(v2xy)):
	if v2flags[i] == 1:
		v2xy[i] = v2xy[i]*(-1.0)	

#convert x,y from arcsec to parsecs:

x = (x / 60 * np.pi / 180)*4600

y = (y / 60 * np.pi / 180)*4600

radius = np.sqrt(x**2 + y**2)

ax2.tricontour(x,y,v2xy, levels=14, linewidths=0.5, colors='k')
cntr2 = ax2.tricontourf(x,y,v2xy, levels=14, cmap="RdBu_r")
fig.colorbar(cntr2,ax=ax2).set_label('v2xy (km/s)')

ax2.set_xlabel('x (arcmin)')
ax2.set_ylabel('y (arcmin)')
fig.suptitle('v2xy')
plt.ylim(-40,40)
plt.xlim(-40,40)
plt.show()

#-----1D radius v dispersions plots


#plt.xlabel('R (pc)')
#plt.ylabel('Velocity Dispersions (km/s)')
#plt.scatter(radius, v2zz, label = "Line of Sight (v2zz)",marker=".")
##plt.scatter(radius, v2xx, label = "PM (v2xx)",marker=".")
##plt.scatter(radius, v2yy, label = "PM (v2yy)",marker=".")
##plt.xscale('log')
#plt.legend()