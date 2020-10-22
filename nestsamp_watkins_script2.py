import numpy as np
# import matplotlib.pyplot as plt
from astropy import table, units as u
import cjam


# import sys

def watkinsMethod(positions, phi, xyFileName, MLratio, beta_num, inclAngle, dataIsRadial, MGEfile,nonLuminousMasses, nonLuminousRadii):
    # read in R values (distance from center of omega-cen) and set up x,y eqns to integrate, units of arcseconds
    Rpos_ac = np.array(positions)

    # convert R from parsecs to radians
    #
    #    Rpos = np.arcsin(Rpos_pc/4600)
    #
    #    #convert R from arcsecs to radians

    Rpos = Rpos_ac * (1 / 3600) * np.pi * (1 / 180)
    x = []
    y = []

    # define x,y relationship to R input and convert from radians to arcseconds

    for i in range(len(Rpos)):
        xpos = Rpos[i] * np.cos(phi) * (180 / np.pi) * 3600
        x.append(xpos)
        ypos = Rpos[i] * np.sin(phi) * (180 / np.pi) * 3600
        y.append(ypos)

    f = open(xyFileName, "w+")
    for i in range(len(Rpos)):
        f.write(str(x[i]) + "," + str(y[i]) + "\n")
    f.close()

    pos = table.QTable.read(xyFileName, format="ascii", names=["x", "y"])
    pos["x"].unit = u.arcsec
    pos["y"].unit = u.arcsec

    # read in the tracer density MGE from the example and add appropriate units
    tracer_mge = table.QTable.read(MGEfile, format="ascii",
                                   names=["n", "i", "s", "q"])
    tracer_mge["i"].unit = 1 / u.pc ** 2
    tracer_mge["s"].unit = u.arcsec

    # set anisotropy for each component of tracer MGE
    beta = np.ones(len(tracer_mge)) * beta_num

    # set rotation for each component of tracer MGE
    kappa = np.array([0, 0.0, -0.4, - 1.1, -0.6, 0, 0, 0])

    # set distance
    distance = 4.6 * u.kpc

    # set inclination
    inclination = inclAngle * u.rad

    # set potential MGE, assumes mass-follows-tracers with mass/tracer ratio of 2.6
    potential_mge = tracer_mge.copy()
    potential_mge["i"] *= MLratio * u.Msun

    # add units to non-luminous components
    nonLuminousMassArray = np.array(nonLuminousMasses) * u.Msun
    nonLuminousRadiusArray = np.array(nonLuminousRadii) * u.arcsec
    
    
    # call CJAM integration to get velocity dispersions
    moments_rot = cjam.axisymmetric(pos["x"], pos["y"], tracer_mge, potential_mge,
                                    distance, nonLuminousMassArray, nonLuminousRadiusArray, beta=beta, kappa=kappa, incl=inclination)

    # define dispersions and velocities
    PMx_Velocities = moments_rot["vx"]
    PMx_VelocityDispersions = moments_rot["v2xx"]
    PMy_Velocities = moments_rot["vy"]
    PMy_VelocityDispersions = moments_rot["v2yy"]

    # convert dispersions from (mas/yr)**2 to (km/s)**2, then find LOS dispersions from the proper motion dispersions
    kms_PMx_VelocityDispersions = (np.sqrt(PMx_VelocityDispersions) * (u.yr) * (1 / (u.mas))) * 21.829
    kms_PMy_VelocityDispersions = (np.sqrt(PMy_VelocityDispersions) * (u.yr) * (1 / (u.mas))) * 21.829
    lineofSightVelocities = moments_rot["vz"]
    lineofSightVelocityDispersions = moments_rot["v2zz"]

    if not dataIsRadial:
        return np.sqrt((1 / 2) * (kms_PMx_VelocityDispersions ** 2 + kms_PMy_VelocityDispersions ** 2))
    else:
        return np.sqrt(lineofSightVelocityDispersions) * u.s * (1/(u.km))
