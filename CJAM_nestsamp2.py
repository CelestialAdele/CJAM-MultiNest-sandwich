#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 15:50:33 2020

@author: addy
"""

import json
import numpy as np
import pymultinest
import nestsamp_watkins_script2 as phi_ws
from scipy.integrate import trapz
from scipy.interpolate import interp1d
import pandas as pd

#see how long this shit takes----------------------------------------------

import atexit
from time import time, strftime, localtime
from datetime import timedelta


def secondsToStr(elapsed=None):
    if elapsed is None:
        return strftime("%Y-%m-%d %H:%M:%S", localtime())
    else:
        return str(timedelta(seconds=elapsed))

def log(s, elapsed=None):
    line = "="*40
    print(line)
    print(secondsToStr(), '-', s)
    if elapsed:
        print("Elapsed time:", elapsed)
    print(line)
    print()

def endlog(start):
    end = time()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))
    

        
# create function that will create list of parameters based on non-luminous matter conditionals   
def createParams(addBlackHole, addDarkMatter, NofDMgaussians):

    paramList = []
    
    if addDarkMatter == True:
        for i in range(NofDMgaussians):
            paramList.append("Log(DM Halo Mass " + str(i+1) + ")")
            paramList.append("Log(DM Halo Radius " + str(i+1) + ")")
        
    return paramList


def GeneratePyMultiNestResults(databaseFileName, outputDirectoryPath, objectName, datasetName, NofLivePoints, Interp, dataIsRadial, MGEfileName, NofInterpPts,addBlackHole,addDarkMatter, NofDMgaussians):
    Interp = Interp
    dataIsRadial = dataIsRadial
    addBlackHole = addBlackHole
    addDarkMatter = addDarkMatter
    
    parameters = createParams(addBlackHole, addDarkMatter, NofDMgaussians)
    
    def prior(cube, ndim, nparams):
    # check if user added DM halo. If true, create a number of gaussian components to describe dm halo equal to 2*NofDMgaussians, all with same prior ranges
        if addDarkMatter == True:
            for i in range(NofDMgaussians):
                # define log dark matter halo mass
                cube[2*i] = cube[2*i] * 4 + 3
                # define log dm halo radius
                cube[2*i+1] = cube[2*i+1] * 2.30 + 1
    
    
    start = time()
    atexit.register(endlog)
    log("Start Program")
    
    # print summary of user inputs
    print("Dataset:" + " " + str(datasetName) + ", Output directory:" + str(outputDirectoryPath) + ", Interpolation:" + str(Interp) + ", N of Interpolated Points:" + str(NofInterpPts) + ", Radial Velocity Dispersion:" + str(dataIsRadial) + ", Black Hole:" + str(addBlackHole) + ", Dark Matter:" + str(addDarkMatter) + ", N of DM gaussians:" + str(NofDMgaussians))
    # store summary of user inputs to file
    with open(outputDirectoryPath + "/summaryOfUserInputs.txt", "w+") as file:
        file.write("Database File Name: " + str(databaseFileName) + ", Dataset:" + " " + str(datasetName) + ", Output directory:" + str(outputDirectoryPath) + ", Interpolation:" + str(Interp) + ", N of Interpolated Points:" + str(NofInterpPts) + ", Radial Velocity Dispersion:" + str(dataIsRadial) + ", Black Hole:" + str(addBlackHole) + ", Dark Matter:" + str(addDarkMatter) + ", N of DM Gaussians:" + str(NofDMgaussians))

    # grab data from globular cluster database
    data = pd.read_pickle(databaseFileName)
    
    # define columns in data

    objectColumn = data["Cluster"]
    RColumn = data["<r>"]
    dispersion_observedColumn = data["σ"]
    up_errorColumn = data["ΔσUp"]
    low_errorColumn = data["ΔσLow"]
    datasetColumn = data["Type"]
    # UNITS: R in arcsecs, dispersion in km/s, dispersion error in km/s

    # make membership cuts:
    R_mem = []
    dispersion_observed_mem = []
    error_mem = []


    for i in range(len(RColumn)):
        if objectColumn[i] == objectName and datasetColumn[i] == datasetName:
            R_mem.append(RColumn[i])
            dispersion_observed_mem.append(dispersion_observedColumn[i])
            error_mem.append((up_errorColumn[i] + low_errorColumn[i])/2)

    # convert lists to numpy array to use in likelihood function
    dispersion = np.array(dispersion_observed_mem)

    # create error array

    error = np.array(error_mem)

    def loglike(cube, ndim, nparams):
        
        
        # attach parameters to prior values, defined in prior fucntion
        
        # append non-luminous parameters to mass and radius lists
        nonLuminousMasses = []
        nonLuminousRadii = []

            
        if addDarkMatter == True:
            if NofDMgaussians == 1:
                dmMass1 = 10**cube[parameters.index("Log(DM Halo Mass 1)")]
                dmRadius1 = 10**cube[parameters.index("Log(DM Halo Radius 1)")]
                nonLuminousMasses.append(dmMass1)
                nonLuminousRadii.append(dmRadius1)
                
            elif NofDMgaussians == 2:
                dmMass1 = 10**cube[parameters.index("Log(DM Halo Mass 1)")]
                dmRadius1 = 10**cube[parameters.index("Log(DM Halo Radius 1)")]
                dmMass2 = 10**cube[parameters.index("Log(DM Halo Mass 2)")]
                dmRadius2 = 10**cube[parameters.index("Log(DM Halo Radius 2)")]
                nonLuminousMasses.append(dmMass1)
                nonLuminousRadii.append(dmRadius1)
                nonLuminousMasses.append(dmMass2)
                nonLuminousRadii.append(dmRadius2)
                
            elif NofDMgaussians == 3:
                dmMass1 = 10**cube[parameters.index("Log(DM Halo Mass 1)")]
                dmRadius1 = 10**cube[parameters.index("Log(DM Halo Radius 1)")]
                dmMass2 = 10**cube[parameters.index("Log(DM Halo Mass 2)")]
                dmRadius2 = 10**cube[parameters.index("Log(DM Halo Radius 2)")]
                dmMass3 = 10**cube[parameters.index("Log(DM Halo Mass 3)")]
                dmRadius3 = 10**cube[parameters.index("Log(DM Halo Radius 3)")]
                nonLuminousMasses.append(dmMass1)
                nonLuminousRadii.append(dmRadius1)
                nonLuminousMasses.append(dmMass2)
                nonLuminousRadii.append(dmRadius2)
                nonLuminousMasses.append(dmMass3)
                nonLuminousRadii.append(dmRadius3)

            
        # get R file values to use as CJAM inputs
        # if Interp = True, interpolate points
        R_values = R_mem
        if Interp == True:
                log_R_values = np.log10(np.array(R_mem))
                evenlyspaced_log_R_values = np.linspace(log_R_values[0], log_R_values[-1], int(NofInterpPts))
                R_values = 10**(evenlyspaced_log_R_values)

        # feed phi_ws a list of phi values:
        phis = np.linspace(0, 2 * np.pi, len(R_values))

        before_averaging_LOSDispersions = []

        # below calls CJAM for len(R_mem) different phi values, which gives back 20 arrays of sigma values, with each array having len(R_mem) R values and 1 phi value. Where len(R_mem) = number of datapoints
        for i in range(len(R_values)):
            LOSdispersions = phi_ws.watkinsMethod(R_values, phis[i], outputDirectoryPath + "/xy.dat", 1.0, 0, 0.8377, dataIsRadial, MGEfileName, nonLuminousMasses, nonLuminousRadii)
            before_averaging_LOSDispersions.append(LOSdispersions)

        # now we need to transpose the sigma arrays to integrate over all phis for one value of R
        transposed_LOSdispersions = []
        for row in range(len(R_values)):
            current_list = []
            for column in range(len(R_values)):
                t_LOS = float(before_averaging_LOSDispersions[column][row])
                current_list.append(t_LOS)
            transposed_LOSdispersions.append(current_list)

        # integrate each transposed_LOS over phi, and then final result is stored in avg_LOSdispersions list
        avg_LOSdispersions = []
        for i in range(len(R_values)):
            averaged_LOSdispersions = (1 / (2 * np.pi)) * trapz(transposed_LOSdispersions[i], x=phis)
            avg_LOSdispersions.append(averaged_LOSdispersions)

        # interpolate then create list of sigmas using the interpolation method
        array_R = np.array(R_values)
        log_R = np.log10(array_R)
        evenlyspaced_log_R = np.linspace(log_R[0],log_R[-1],len(R_values))
        f = interp1d(evenlyspaced_log_R, avg_LOSdispersions, kind="cubic")

        # find new sigma values based on interpolation results and store in dispersion_CJAM
        array_R_mem = np.array(R_mem)
        log_R_mem = np.log10(array_R_mem)
        dispersion_CJAM = []
        for i in range(len(log_R_mem)):
            newsigma = f(log_R_mem[i])
            dispersion_CJAM.append(newsigma)

        # define likelihood function from eqn 29 in the following paper: https://arxiv.org/pdf/1805.05883.pdf
        likelihood = (2 * np.pi * error ** 2) ** (-1 / 2) * np.exp(-(dispersion - dispersion_CJAM) ** 2 / (2 * error ** 2))

        return np.log(likelihood).sum()


    ##------------------------------------------------------------------------------
        
    # number of dimensions our problem has
    
    n_params = len(parameters)

    # run MultiNest
    pymultinest.run(loglike, prior, n_params, outputfiles_basename = outputDirectoryPath + "/",
        resume = False, verbose = True, n_live_points=NofLivePoints)

    json.dump(parameters, open(outputDirectoryPath + '/params.json', 'w')) # save parameter names

     #now run the script and analyse the output using multinest_marginals.py::

      #  $ python 1_1d_multimodal.py &&  multinest_marginals.py 1_1d_multimodal_out

     #then open the file 1_1d_multimodal_outmarg.pdf

     #Btw, ln(ev) should be ln(1 / 2)


    #--------------------------------------------------------------------------------------
     #see how long this shit takes
    with open(outputDirectoryPath + "/summaryOfUserInputs.txt", "a") as file:
        file.write('\n' + str(endlog(start)))
    endlog(start)

