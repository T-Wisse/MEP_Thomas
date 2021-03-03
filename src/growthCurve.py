# -*- coding: utf-8 -*-
"""
Processing lab data:
    Using spectrophotometer to on different dilutions and timepoints to create a growth curve
    Fit a line to the growth curve to determine slope (=doubling time)   

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
os.chdir('D:/Users/Thomas/Studie/MEP/MEP_Thomas/src')

#%% 01/02/21 YLL3a in YPD Dilution

OD = [2.410, 1.831, 0.955, 0.542, 0.276]
dilution = [1, 2, 5, 10, 20]

plot = plt.scatter(dilution, OD)
plt.title('OD of dilution of YLL3a in YPD')
plt.xlabel('Dilution factor')
plt.ylabel('OD')
            

#%% 11/02/21 growth curves of YLL3a in YPD and CSM



def func(x, a, b):
     return 2**(a*x+b)

dataYPD = pd.read_excel('../data/growthCurveYPD210211.xlsx')
dataCSM = pd.read_excel('../data/growthCurveCSM210211.xlsx')

tYPD = np.array(np.cumsum(dataYPD['Time']), dtype = float)
odYPD = np.array(dataYPD['YPD OD']) * np.array(dataYPD['dilution'], dtype = float)

tYPDlinear = tYPD[tYPD>65] #linear region in YPD starts around 50 min
odYPDlinear = odYPD[tYPD>65]

tCSM = np.array(np.cumsum(dataCSM['Time']), dtype = float)
odCSM = np.array(dataCSM['CSM OD'])

tCSMlinear = tCSM[tCSM>100] #linear region in CSM starts around 100 min
odCSMlinear = odCSM[tCSM>100] 

varsYPD, pcov = curve_fit(func, tYPDlinear, odYPDlinear, p0 = [0.005,1], bounds  = (-10, 10))
varsCSM, pcov = curve_fit(func, tCSMlinear, odCSMlinear, p0 = [0.005,1], bounds  = (-10, 10))

dtYPD = 1/varsYPD[0] # doubling time YPD
dtCSM = 1/varsCSM[0] # doubling time CSM


print("Doubling time for WT Yeast in YPD is " + str(int(round(dtYPD))) + " minutes")
print("Doubling time for WT Yeast in CSM is " + str(int(round(dtCSM))) + " minutes")

#%% calculate R squared
residuals = odYPD - func(tYPD, *varsYPD)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((odYPD-np.mean(odYPD))**2)
r_squared = 1 - (ss_res / ss_tot)

#%% plot growth curvep
plt.plot(tYPD,odYPD, label = 'YPD data')
plt.plot(tCSM,odCSM, label = 'CSM data')
plt.xlabel('Time incubated (min)')
plt.ylabel('OD')
plt.title('Growth curve of WT yeast in YPD and CSM')

# labelYPDfit = 'YPD fit, doubling time:' + str(int(round(1/dtYPD))) + ' min'
# labelCSMfit = 'CSM fit, doubling time:' + str(int(round(1/dtCSM))) + ' min'

labelYPDfit = 'YPD fit'
labelCSMfit = 'CSM fit'

plt.plot(tYPD, func(tYPD, *varsYPD), label = labelYPDfit)
plt.plot(tCSM, func(tCSM, *varsCSM), label = labelCSMfit)
plt.legend()
plt.yscale('log', basey = 2)
plt.grid()