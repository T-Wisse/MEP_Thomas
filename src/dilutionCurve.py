# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 16:22:14 2021

@author: Thomas
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
os.chdir('D:/Users/Thomas/Studie/MEP/MEP_Thomas/src')
#%%

dataOD = pd.read_excel('../../data/labdata/dilutionCurve110221.xlsx')

plt.plot(dataOD['Dilution factor'], dataOD['CSM OD'], label = 'OD CSM')
plt.plot(dataOD['Dilution factor'], dataOD['YPD OD'], label = 'OD YPD')
plt.legend()
plt.title('OD of dilution curve yeast in CSM and YPD')
plt.xlabel('Dilution factor')
plt.ylabel('OD')
plt.xscale('log', basex = 2)
plt.yscale('log', basey = 2)
plt.grid()
