#The purpose of this code is to be able to use data from MW globular clusters with only velocity dispersions and distances from the center of the cluster. This code sets up an integration of the LOS dispersions outputted by CJAM over all angles, so that we get averaged LOS dispersions from "R values" (distance from the center). We also interpolate additional dispersions in this script using scipy interpolation functions. - 5/26/2020

#kind = linear, cubic

import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import phi_integral_watkins_script as phi_ws
import random

#feed phi_ws a list of phi:
R = np.loadtxt("Rdat20.dat")

phis = np.linspace(0,2*np.pi,len(R))



LOSdispersions = []
trialLOS = []

#below calls CJAM for 20 different phi values, which gives back 20 arrays of sigma values, with each array having 20 R values and 1 phi value
for i in range(len(R)):
    LOSdispersions = phi_ws.watkinsMethod(R,100000,1,phis[i])
    trialLOS.append(LOSdispersions)


#now we need to transpose the sigma arrays to integrate over all phis for one value of R
transposed_LOS = []
for row in range(len(R)):
    current_list = []
    for column in range(len(R)):
        t_LOS = float(trialLOS[column][row])
        current_list.append(t_LOS)
    transposed_LOS.append(current_list)
    
#integrate each transposed_LOS over phi
avg_LOSdispersions = []
for i in range(len(R)):
    averaged_LOSdispersions = (1/(2*np.pi))*trapz(transposed_LOS[i], x = phis)
    avg_LOSdispersions.append(averaged_LOSdispersions)

#-----------------------------------------------------------------

#interpolate this ish
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
f = interp1d(R, avg_LOSdispersions, kind = "cubic")
newR = np.arange(R[0],R[19],1)
new_sigma = f(newR)

#-----------------------------------------------------------------

#take 10 R points from interpolation and get new dispersions, plot against calculated dispersions

ten_Rvalues = random.choices(newR, k = 10)
ten_Rvalues.sort()

#ten_Rvalues = [np.log10(1),np.log10(2),np.log10(3),np.log10(4),np.log10(5),np.log10(6),np.log10(7),np.log10(8),np.log10(9),np.log10(10)]
ten_dispersions = f(ten_Rvalues)

#run CJAM with R values from randomly selected interpolated points. Then we can better compare the interpolation with CJAM calculations.

#change back from logR to R so that CJAM can correctly do its thang
#a_ten_Rvalues = np.array(ten_Rvalues)
#a_ten_Rvalues = 10**(a_ten_Rvalues)

f = open("newRvalues.dat", "w+")
for i in range(len(ten_Rvalues)):
    f.write(str(ten_Rvalues[i]) + "\n")
f.close()

#feed phi_ws a list of phi:
checkR = np.loadtxt("newRvalues.dat")
checkphis = np.linspace(0,2*np.pi,len(checkR))

cLOSdispersions = []
ctrialLOS = []

#below calls CJAM for 20 different phi values, which gives back 20 arrays of sigma values, with each array having 20 R values and 1 phi value
for i in range(len(checkR)):
    cLOSdispersions = phi_ws.watkinsMethod(checkR,100000,1,checkphis[i])
    ctrialLOS.append(cLOSdispersions)


#now we need to transpose the sigma arrays to integrate over all phis for one value of R
ctransposed_LOS = []
for row in range(len(checkR)):
    ccurrent_list = []
    for column in range(len(checkR)):
        ct_LOS = float(ctrialLOS[column][row])
        ccurrent_list.append(ct_LOS)
    ctransposed_LOS.append(ccurrent_list)
    
#integrate each transposed_LOS over phi
cavg_LOSdispersions = []
for i in range(len(checkR)):
    caveraged_LOSdispersions = (1/(2*np.pi))*trapz(ctransposed_LOS[i], x = checkphis)
    cavg_LOSdispersions.append(caveraged_LOSdispersions)


#a_checkR = np.array(checkR)
#
#a_checkR = 10**(a_checkR)
#-----------------------------------------------------------------

#plots n stuff    

#plt.scatter(R, avg_LOSdispersions, color = 'darkorchid', label = "CJAM Averaged LOS Dispersions (km/s)" )
#plt.scatter(ten_Rvalues,ten_dispersions, color = 'darkorange', label = "Interpolated LOS Dispersions (km/s)")
plt.scatter(ten_Rvalues, ten_dispersions, label = "initial interpolation results", color = 'darkorange')
plt.scatter(checkR, cavg_LOSdispersions, label = "calculated CJAM results", color = 'darkorchid')
plt.legend()
plt.xlim(0,50)
plt.show()
