#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:36:18 2020

@author: andrewmartin

fitShrinkSwell.py

"""

import shrinkSwellTools as sst
#import matplotlib.pyplot as plt
import numpy as np


#
# SET THE TYPE OF SIMULATION
#
# 1 - shrink only simulation (e.g. PBS)
#
# 2 - shrink & swell  (i.e. with solute)
#
shrinkswell = 2

#
# Fit a curve or just simulate it?
#
# 1- Simulate a curve (single set of parameters)
# 2- Fit a parameter to the data (which parameter is fitted depends on shrinkswell) 
#
simORfit = 2

#
# create an object with all the fitting routines & parameters
#
sstb = sst.shrinkSwellToolBox(ss=shrinkswell)

#
# Choose whether to use data (True or False)
#
usedata = True

#
# set filename of data to read if required
#
# datapath should end with a slash / or \ depending on the operating system
#
# (data type is assumed to be a .csv file with two columns: time & volume)
# time in seconds; Volume units are arbitrary (normalised to starting volume)
#
# (set datapath to None for current working directory)
#
if usedata == True:
    datapath = "/Users/Saffron/Dropbox/Saffrons stuff/Research/RMIT/Experimental/PythonFitting/1pt11_DES/"
    data_filename = "ChClGly2.csv"
    

#
# Set parameters for output files
#
# prefix will begin filenames for all output files (use it not to overwrite results)
#
sstb.outpath = "/Users/Saffron/Dropbox/Saffrons stuff/Research/RMIT/Experimental/PythonFitting/1pt11_DES/"
sstb.prefix = "ChClGly2_Vb.2789_Lp401_PSVary"
plot_ext = ".png"  #type of file for output images


#
# extracellular osmolality (or final osmolality)
#
# total:
sstb.Ce = 1663*1e-3 

# solute:
sstb.CeS = 1368*1e-3 

# Starting osmolality
sstb.C0 = 295*1e-3  # mosm/kg 

# initial internal solute osmolality
sstb.CiS = 0.0


#
# Volume parameters
#
# cell radius
r = 10.0 # cell radius in micron


# Vb as a fraction of V0 (adjusted below)

sstb.Vb= 0.2789

# partial molar fraction
sstb.VbarS =  .7103 #0.85943 # L/mol

# initial solute volume
sstb.Vs = 0.0 


#
# Lp controls the shrinking rate
# - this value is used for simulated plots and fitting the swell parameter (Ps)
#
Lpfixed = .401/60        # micron/min/atm  converted to seconds 

#
# Ps 
#- this value is only use for simulated plots (not fitting)
#
Psfixed = 0.001/60       # micron / min   /60 converted to seconds

#
# Fitting parameters for Lp  (shrink parameter)
# 
# - only used if shrinkswell = 1 AND simORfit = 2
#
# - Start and End define the range of values to test
# -  N is the number of sampling points in the range
#
sstb.LpStart = .001 / 60.
sstb.LpEnd = 100 /60.
sstb.LpN = 1000

#
# Fitting parameters for Ps   (swell parameter)
# 
# - only used if shrinkswell = 2 AND simORfit = 2
#
# - Start and End define the range of values to test
# -  N is the number of sampling points in the range
#
sstb.PsStart = 0.005 / 60.
sstb.PsEnd = 10 / 60.
sstb.PsN = 1000

#
# Set the time samnpling
#
# If total_time exceeds the max value in the data file it will be overwritten
#
# nt is the number of time points
# 
sstb.total_time = 1200
sstb.nt=600


# Temperature
sstb.T = 293   #Kelvin




#
# SOME DERIVED VALUES & ADJUSTMENTS
#

# cell membrane area (micron^2)
sstb.A =  4 * np.pi * (r**2)    #10 micron radius
# initial cell volume:
sstb.V0 = 1.0 * (4/3)*np.pi*(r**3)

# initial value of Volume Vc
sstb.Vc = 1.0 * sstb.V0

sstb.Vb *= sstb.V0


sstb.shrinkswell = shrinkswell
sstb.simORfit = simORfit
sstb.usedata = usedata
if usedata==True:
    sstb.datapath = datapath
    sstb.data_filename = data_filename
sstb.plot_ext = plot_ext

sstb.Lp = Lpfixed
sstb.Ps = Psfixed
sstb.radius = r

#
# generate the tag
# 
sstb.tag = sstb.prefix
if usedata==True:
    sstb.tag += "_"+data_filename[:-4]
    
if shrinkswell == 1:
    sstb.tag += "_shrinkonly"
elif shrinkswell == 2:
    sstb.tag += "_shrinkswell"
    
if simORfit == 1:
    sstb.tag += "_sim"
elif simORfit == 2:
    sstb.tag += "_fit"


#################################################
#                                               #
#     All parameters have been set              #
#                                               #
#     Run the analysis!                         #
#                                               #
#################################################

#
# write out all the parameters to file
#
sstb.write_all_params_to_file()



sstb.run() #shrinkswell,simORfit,usedata,datapath,data_filename,plot_ext)
