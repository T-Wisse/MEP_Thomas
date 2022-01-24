# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 00:52:23 2021

@author: Thomas

plotting lab results
"""

import matplotlib.pyplot as plt

OD = [2.410, 1.831, 0.955, 0.542, 0.276]
dilution = [1, 2, 5, 10, 20]

plot = plt.scatter(dilution, OD)
plt.title('OD of dilution of YLL3a in YPD')
plt.xlabel('Dilution factor')
plt.ylabel('OD')
            


