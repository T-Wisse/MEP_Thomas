# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:57:47 2020

@author: Thomas
"""

from timeit import default_timer as timer

import numpy as np
from collections import defaultdict 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os
os.chdir('D:/Users/Thomas/Studie/MEP/MEP_Thomas/src')

path = 'D:/Users/Thomas/Studie/MEP/MEP_Thomas/src/python_modules'
import sys
sys.path.insert(0, path) #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

from modulesT import common_interactors_T
from modulesT import common_go_children
from modulesT import getChildrenGoTerm

# from python_modules.module_common_measures import common_partners

#%% getting the data

data_go = pd.read_excel('../data/slim-goterms-filtered-data.xlsx')
data_go = data_go.dropna()

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type' ]

data=pd.read_excel('../data/data-BioGrid-Yeast.xlsx')
#%% query
query=['BEM2'] #should be 381 for BEM1 ???


#%% Calling the function common_interactors
# start = timer()
commonInteractorData=common_interactors_T(query,data)
# end = timer()
# print(end-start)

#common_go=common_go(data_go=data_go,data_common_partners=common_partners_data)

#%% get all unqiue interactor genes
uniqueInteractorGenes = list(set(commonInteractorData['interactorGene']))
# uniqueInteractorGenes = commonInteractorData.interactorGene.unique()
#%% for each unqiue interactor gene get its go terms (dict)

goTermsInteractors = []
for ii in uniqueInteractorGenes:
    # print(ii)
    goTermsInteractors.append(list(data_go['go-term'][data_go['Gene'].str.match(ii)]))


#%% for each go term of each interactor gene, get its 1st layer children
# for each interactor, find the overlapping set of 1st layer GO terms with our query gene

print('using cell polarity as parent query GO term')
childrenGoTermsQuery = getChildrenGoTerm('establishment or maintenance of cell polarity') #should we do this for every query gene? check plan on github
commonChildGoTerms = defaultdict(dict)
jj = 0

for ii in uniqueInteractorGenes: #this loop takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?
    print(jj)
    queryGoTerm = goTermsInteractors[jj] #discuss with Leila what to do here...
    for kk in queryGoTerm:
        childrenGoTermsInteractor = getChildrenGoTerm(kk) 
        tempCommonChildGoTerms = childrenGoTermsQuery + childrenGoTermsInteractor
        commonChildGoTerms[ii] = set([x for x in tempCommonChildGoTerms if tempCommonChildGoTerms.count(x) > 1]) #finds all common go terms - empty so far, does it work properly?
    jj = jj+1
# GOTermsInteractors = getGOTermsGene() 

#%%


#%% testing
# childrenGOtermsInteractors = getChildrenGoTerm('establishment or maintenance of cell polarity')


# commonInteractors = defaultdict(dict)
# geneList = ['BEM1', 'BEM3']
# query = ['BEM2', 'BEM3']

# for i in range(1):
#     for genes in geneList:
#         commonInteractors[genes]['query']=query[i]
#         commonInteractors[genes]['names of genes']=genes


# interactionData = data
# queryGene = 'BEM1'

# query1 = interactionData[interactionData['gene-query-name']==queryGene] #Reduce list to only query gene
# interactorList = query1['gene-target-name'].unique() #List genes interacting with query1 gene 

# query1['gene-target-name'] in interactorList


# query1[pd.DataFrame(query1.'gene-target-name'.tolist()).isin(interactorList).any(1).values]

# queryGene = 'a'
# interactorGene = 'b'
# tempCommonInt = 'c'


# commonInteractors[queryGene+'-'+interactorGene] = [queryGene]
# commonInteractors[queryGene+'-'+interactorGene] = [interactorGene]
# commonInteractors[queryGene+'-'+interactorGene] = [tempCommonInt]


