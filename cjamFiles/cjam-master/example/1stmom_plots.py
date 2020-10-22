import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



vx,vy,vz = np.loadtxt('nc_new_moments.dat', dtype='float', usecols=(0,1,2), unpack = True)

x,y =  np.loadtxt('xy.dat', dtype='float',delimiter=",", usecols=(0,1), unpack = True)

#fig =plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#bx = fig.add_subplot(111, projection='3d')
#cx = fig.add_subplot(111, projection='3d')

#ax.scatter(x,y,vx, marker='o')
#plt.show()

#bx.scatter(x,y,vy, marker='o')
#cx.scatter(x,y,vz, marker='o')

#plt.show()

x = x*0.0167
y = y*0.0167

fig, ax2 = plt.subplots()

ax2.plot(x,y, 'ko', ms=3, alpha = 0.1)





#vy =vy*4.74*5.0

		
vx = vx*4.74*5.0
vy = vy*4.74*5.0

ax2.tricontour(x,y,vz, levels=14, linewidths=0.5, colors='k')
cntr2 = ax2.tricontourf(x,y,vz, levels=14, cmap="RdBu_r")
fig.colorbar(cntr2,ax=ax2).set_label('vz (km/s)')
plt.ylim(-20,20)
plt.xlim(-20,20)
plt.title("vz")
ax2.set_xlabel('x (arcmin)')
ax2.set_ylabel('y (arcmin)')
plt.show()
