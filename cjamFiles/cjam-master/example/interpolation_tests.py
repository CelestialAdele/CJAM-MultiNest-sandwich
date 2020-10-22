#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:04:25 2020

@author: addy
"""

import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import phi_integral_watkins_script as phi_ws


#feed phi_ws a list of phi:
R = np.loadtxt("Rdat.dat")
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

plt.scatter(R, avg_LOSdispersions)
plt.ylim(bottom=6)
