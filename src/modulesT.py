# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:58:56 2020

@author: Thomas
"""

import numpy as np
import pandas as pd
# from collections import defaultdict 

#%% Common interactors: find all genes known to have an interaction with both the query gene and its interactors
def common_partners_T(query,interactionData): #interaction data must have a 'gene query name' and a 'gene target name' column
    # d2 = defaultdict(dict)
    
    # testQ = interactionData['gene-query-name']
    
    for geneName in query: #For each query gene
    
        q1 = interactionData[interactionData['gene-query-name']==geneName] #Reduce list to only query gene
        interactorList=q1['gene-target-name'].unique() #List genes interacting with query gene 
        
        for interactor in interactorList:
            
            q2 = interactionData[interactionData['gene-query-name']==interactor]
            interactors=q2['gene-target-name'].unique()
                
            commonInteractors = []
            for interactorName  in interactors:
                if interactorName in interactorList: # if a gene interactor of the query1 is in interactors of query 2 
                    commonInteractors.append(interactorName) #add gene to list
                    
    return commonInteractors #note, so far will only work properly for 1 gene in query...

def gaussFit:
   

# def common_go(data_go,data_common_partners):

                    