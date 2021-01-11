# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:57:47 2020

@author: Thomas
"""

import numpy as np
from collections import defaultdict 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir('D:/Users/Thomas/Studie/MEP/MEP_Thomas/src')
from modulesT import common_partners_T

# from python_modules.module_common_measures import common_partners

#%% getting the data

data_go=pd.read_excel('../datasets/slim-goterms-filtered-data.xlsx')

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type' ]

data=pd.read_excel('../datasets/data-BioGrid-Yeast.xlsx')
#%% query
#query=['BEM2']
query=['ACT1']

#%% Calling the function common_partners

common_partners_data=common_partners_T(query,data)

#common_go=common_go(data_go=data_go,data_common_partners=common_partners_data)