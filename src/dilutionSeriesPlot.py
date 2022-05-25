# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 00:52:23 2021

@author: Thomas

plotting lab results
"""

#%% Old dilution series plot

# import matplotlib.pyplot as plt

# OD = [2.410, 1.831, 0.955, 0.542, 0.276]
# dilution = [1, 2, 5, 10, 20]

# plot = plt.scatter(dilution, OD)
# plt.title('OD of dilution of YLL3a in YPD')
# plt.xlabel('Dilution factor')
# plt.ylabel('OD')
            
#%%

import matplotlib.pyplot as plt
import pandas as pd

path = 'D:/Users/Thomas/Studie/MEP/Data/LabData/'
filename = 'dilutionCurve110221.xlsx'

df = pd.read_excel(path+filename)

dilutionFactor = df['Dilution factor'].tolist()
csmOD = df['CSM OD'].tolist()
ypdOD = df['YPD OD'].tolist()

plt.rc('font', size=16)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
fig, ax = plt.subplots(dpi = 200)



ax.semilogx(dilutionFactor, ypdOD, label = 'YPD')
ax.semilogx(dilutionFactor, csmOD, label = 'CSM')
ax.legend()


plt.title('OD of dilution series of YLL3a in different media')
plt.xlabel('Dilution factor')
plt.ylabel('OD')
