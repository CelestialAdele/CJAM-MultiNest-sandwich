#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:14:30 2020

@author: addy
"""

f = open("Rdat25.dat", "w+")
for i in range(25):
    f.write(str(i*2) + '\n')
f.close()