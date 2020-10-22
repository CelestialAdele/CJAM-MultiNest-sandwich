#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
from astropy import table, units as u
import cjam
#import sys

def watkinsMethod(positions, BHmass, BHradius, phi):
    
    #read in R values (distance from center of omega-cen) and set up x,y eqns to integrate
    
    Rpos_pc = positions
    #convert R from pc to radians
    Rpos = np.arcsin(Rpos_pc/4600)
    x = []
    y = []
    
    #define x,y relationship to R input and convert from radians to arcseconds
    
    for i in range(len(Rpos)):
        xpos = Rpos[i]*np.cos(phi)*(180/np.pi)*3600
        x.append(xpos)
        ypos = Rpos[i]*np.sin(phi)*(180/np.pi)*3600
        y.append(ypos)
    
    f = open("phi_xy.dat", "w+")
    for i in range(len(Rpos)):
        f.write(str(x[i]) + "," + str(y[i]) + "\n")
    f.close()
    
    pos = table.QTable.read("phi_xy.dat", format = "ascii", names = ["x","y"])
    pos["x"].unit = u.arcsec
    pos["y"].unit = u.arcsec
    
    # read in the tracer density MGE from the example and add appropriate units
    tracer_mge = table.QTable.read("paper_mge.dat", format="ascii",
        names=["n","i","s","q"])
    tracer_mge["i"].unit = 1/u.pc**2
    tracer_mge["s"].unit = u.arcsec
    
    # set anisotropy for each component of tracer MGE
    beta = np.ones(len(tracer_mge))*0.01
    
    # set rotation for each component of tracer MGE
    kappa = np.array([0, 0.0, -0.4,- 1.1, -0.6, 0, 0, 0])
    
    # set distance
    distance = 4.6*u.kpc
    
    # set inclination
    inclination = 0.87*u.rad
    
    # set potential MGE, assumes mass-follows-tracers with mass/tracer ratio of 2.6
    potential_mge = tracer_mge.copy()
    potential_mge["i"] *= 2.6*u.Msun
    

    
    #write BHM and BHR as command line arguments
    BH_mass = int(BHmass)*u.Msun
    BH_radius = int(BHradius)*u.arcsec
    
    #call CJAM integration to get velocity dispersions
    moments_rot = cjam.axisymmetric(pos["x"], pos["y"], tracer_mge, potential_mge,
        distance, beta=beta, kappa=kappa, incl=inclination, mbh=BH_mass,
        rbh=BH_radius)
    
    #define dispersions and velocities
    lineofSightVelocities = moments_rot["vz"]
    lineofSightVelocityDispersions = moments_rot["v2zz"]

    # get v2zz and get rid of rid of units so it can be used for pymultinest nested sampling
    return np.sqrt(lineofSightVelocityDispersions) * u.s * (1/(u.km))



