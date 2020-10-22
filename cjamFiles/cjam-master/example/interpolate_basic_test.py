#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:57:11 2020

@author: addy
"""
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#test scipy interpolate function with sine wave

#interpolate this ish

x = np.linspace(0,2*np.pi,5)

y = np.sin(x)

f = interp1d(x,y)
newx = np.linspace(0.1,2*np.pi,10)
new_y = f(newx)



test_y = np.sin(newx)

plt.scatter(newx, new_y, label = "interpoaltion results")
plt.scatter(newx, test_y, label = "calculation results")
plt.legend()
plt.show()