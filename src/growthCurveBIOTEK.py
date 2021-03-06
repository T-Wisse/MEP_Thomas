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
import datetime
 #change to dir of folder

fileDir = os.path.dirname(os.path.abspath(__file__)) # get path to this file
os.chdir(fileDir) # Change working directory to current file
dataDir = '../data/' #path to data relative from filelocation

#%% 18/04/21 load data

datafile = dataDir + 'BIOTEK_210415_yTW001ABC_and_controls_trimmed.xlsx'#Select Biotek data file
dataBiotek = pd.read_excel(datafile) #load Biotek data
# dataBiotek = dataBiotek.dropna() #If datafile is not trimmed

# titles: row 28
# data: rows 29 - 366
        
#%% For each strain get all data of a well as np array

rows = ['B','C','D','E','F','G']

totalTime3a = 0
totalTime117 = 0
totalTimeTWA = 0
totalTimeTWB = 0
totalTimeTWC = 0

for row in rows: #Go through the different rows (B through G)
    yTWA = np.array(dataBiotek[row +'2']) #Select column
    yTWB = np.array(dataBiotek[row +'4'])
    yTWC = np.array(dataBiotek[row +'6'])
    y3a = np.array(dataBiotek[row +'8'])
    y117 = np.array(dataBiotek[row +'10'])
    yNegControl = np.array(dataBiotek[row +'12'])
    
    
    #%% get time array
    t = np.array(dataBiotek['Time']) #raw time data from Biotek
    
    time = [None] * len(t) #initialize time list
    
    # CAN WE DO THIS WITHOUT A FOR LOOP? SHOULD BE POSSIBLE FOR NP.ARRAY ...
    for ii in list(range(0,171)): #Change time format to seconds (1/2)
        time[ii] = (datetime.timedelta(hours=t[ii].hour, minutes=t[ii].minute, seconds=t[ii].second).total_seconds()) #insane workaround because values below 1 day do not have the day attribute....
        
    for ii in list(range(171,len(t))): #Change time format to seconds (2/2)
        time[ii] = (datetime.timedelta(days=t[ii].day, hours=t[ii].hour, minutes=t[ii].minute, seconds=t[ii].second).total_seconds())
            
    time = np.array(time)/60 #convert time to minutes
    
    # t_temp = np.array(range(0,338))
    # t_approx = ((t_temp*8.245)+3.34) #t in minutes. Approximation of acutal time... Should look into how I get time in seconds from a datatime.time format
            
    #%% get linear regime data for determining doubling time and fit func. Use fit to find doubling time
    #linear regime determined manually
    
    #Fitting function
    def func(x, a, b): 
         return 2**(a*x+b)
     
    #get only linear data
    y3aLinear = y3a[45:66] 
    y117Linear = y117[52:73]
    yTWALinear = yTWA[73:95]
    yTWBLinear = yTWB[73:95]
    yTWCLinear = yTWC[95:123]
    
    #get time corresponding to linear data
    t3aLinear = time[45:66] 
    t117Linear = time[52:73]
    tTWALinear = time[73:95]
    tTWBLinear = time[73:95]
    tTWCLinear = time[95:123]
    
    #perform fit using func and linear data
    vars3a, pcov = curve_fit(func, t3aLinear, y3aLinear, p0 = [0.005,1], bounds  = (-10, 10)) 
    vars117, pcov = curve_fit(func, t117Linear, y117Linear, p0 = [0.005,1], bounds  = (-10, 10))
    varsTWA, pcov = curve_fit(func, tTWALinear, yTWALinear, p0 = [0.005,1], bounds  = (-10, 10))
    varsTWB, pcov = curve_fit(func, tTWBLinear, yTWBLinear, p0 = [0.005,1], bounds  = (-10, 10))
    varsTWC, pcov = curve_fit(func, tTWCLinear, yTWCLinear, p0 = [0.001,1], bounds  = (-10, 10))
    
    #get doubling time from fit
    dt3a = 1/vars3a[0]
    dt117 = 1/vars117[0]
    dtTWA = 1/varsTWA[0]
    dtTWB = 1/varsTWB[0]
    dtTWC = 1/varsTWC[0]
    
    print("Doubling time for yLL3a is " + str(int(round(dt3a))) + " minutes")
    print("Doubling time for yLL117 is " + str(int(round(dt117))) + " minutes")
    print("Doubling time for yTW001A is " + str(int(round(dtTWA))) + " minutes")
    print("Doubling time for yTW001B is " + str(int(round(dtTWB))) + " minutes")
    print("Doubling time for yTW001C is " + str(int(round(dtTWC))) + " minutes")
    
    #Save the total time for later averaging 
    totalTime3a = totalTime3a + int(round(dt3a))
    totalTime117 = totalTime117 + int(round(dt117))
    totalTimeTWA = totalTimeTWA + int(round(dtTWA))
    totalTimeTWB = totalTimeTWB + int(round(dtTWB))
    totalTimeTWC = totalTimeTWC + int(round(dtTWC))
    
    
    #%% plot data
    
    
    plt.figure(dpi=500) #set Figure quality
    plt.yscale('log', basey = 2)
    
    #plot Biotek data
    plt.plot(time, yTWA, label = 'yTWA')
    plt.plot(time, yTWB, label = 'yTWB')
    plt.plot(time, yTWC, label = 'yTWC')
    plt.plot(time, y3a, label = 'yLL3a')
    plt.plot(time, y117, label = 'yLL117')
    plt.plot(time, yNegControl, label = 'Neg Control')
    
    #Plot fitting result
    plt.plot(t3aLinear, func(t3aLinear, *vars3a), label = '3a Fit')
    plt.plot(t117Linear, func(t117Linear, *vars117), label = 'yLL117 Fit')
    plt.plot(tTWALinear, func(tTWALinear, *varsTWA), label = 'yTW001A Fit')
    plt.plot(tTWBLinear, func(tTWBLinear, *varsTWB), label = 'yTW001B Fit')
    plt.plot(tTWCLinear, func(tTWCLinear, *varsTWC), label = 'yTW001C Fit')
    
    plt.xlabel('Time incubated (minutes)')
    plt.ylabel('OD')
    plt.title('Row' + row + 'Growth curves of yTW001 and controls in YPD')
    
    plt.legend()
    plt.grid()


# labelYPDfit = 'YPD fit, doubling time:' + str(int(round(1/dtYPD))) + ' min'
# labelCSMfit = 'CSM fit, doubling time:' + str(int(round(1/dtCSM))) + ' min'

# labelYPDfit = 'YPD fit'
# labelCSMfit = 'CSM fit'


# plt.plot(tCSM, func(tCSM, *varsCSM), label = labelCSMfit)

#%%


print("Average doubling time for yLL3a is " + str(totalTime3a/6) + " minutes")
print("Average doubling time for yLL117 is " + str(totalTime117/6) + " minutes")
print("Average doubling time for yTW001A is " + str(totalTimeTWA/6) + " minutes")
print("Average doubling time for yTW001B is " + str(totalTimeTWB/6) + " minutes")
print("Average doubling time for yTW001C is " + str(totalTimeTWC/6) + " minutes")