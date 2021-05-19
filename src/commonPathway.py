# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:52:03 2021

@author: Thomas

Goal: get relation between pathways of genes

"""


import pandas as pd
import numpy as np
import os as os
os.chdir('D:/Users/Thomas/Studie/MEP/MEP_Thomas/src')
from python_modules.modulesT import getPathwayForGeneFunc as gp

#%% get pathway data on query
os.chdir('D:/Users/Thomas/Studie/MEP/python-modules-for-bioinformatic-analyses/src')
query = pd.read_excel('../datasets/genesCellPolarity_SGD_amigo.xlsx')
query = np.unique(query).tolist()

# query = ['fas1','fas2']            
a = gp(query)

#%%
