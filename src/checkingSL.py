# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:07:05 2021

@author: Thomas

get common interactors for query and its synthetic lethal interactors
"""
import pandas as pd
import numpy as np

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return(a_set.intersection(b_set))  
    return(False)

query = pd.read_excel('../datasets/genesCellPolarity_SGD_amigo.xlsx')
query = np.unique(query).tolist()
bem1SLint=pd.read_excel('../datasets/bem1interactorsSL.xlsx')

bem1SLint = bem1SLint['Gene'].tolist()
print(common_member(bem1SLint, query))