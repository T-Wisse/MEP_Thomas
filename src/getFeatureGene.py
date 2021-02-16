# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 10:36:50 2021

@author: Thomas

Central script to call functions which import data from yeastmine and save to excel

"""
import pathlib
loc = pathlib.Path(__file__).parent.absolute() #get folder location
import os
os.chdir(loc) #move to correct folder...
import numpy as np
import pandas as pd
import python_modules.modulesT as mt

#%% set tile of file
title = 'myYeastDataSet'


#%% getting list query genes

query = pd.read_excel('../datasets/genesCellPolarity_SGD_amigo.xlsx')
query = np.unique(query).tolist()

#%% call function from intermine

# query = 'BEM1'
feature = mt.getPathwayForGene(query)

#%% save
fileDir = os.path.dirname(os.path.realpath('__file__'))
try:
    fileLoc = os.path.join(fileDir, '../datasets')
except:
    print('Could not find dataset folder, saving in current folder')
    fileLoc = fileDir
finally:
    df = pd.DataFrame(feature)
    writer = pd.ExcelWriter(title + '.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='yeastmine', index=False)
    writer.save()
    