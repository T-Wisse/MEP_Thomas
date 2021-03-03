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
query=['BEM1', 'BEM2'] #should be 381 for BEM1 ???


#%% Calling the function common_interactors
# start = timer()
commonInteractorData=common_interactors_T(query,data)
# end = timer()
# print(end-start)

#common_go=common_go(data_go=data_go,data_common_partners=common_partners_data)

#%% get all unqiue interactor genes
uniqueInteractorGenes = list(set(commonInteractorData['interactorGene']))
# uniqueInteractorGenes = commonInteractorData.interactorGene.unique()
#%% for each unqiue interactor gene get its go terms

def goTermsFromGenes(genes,data_go):
    '''
    Returns all go terms from one or more genes based on supplied data_go
    '''
    
    goTerms = []
    for ii in genes:
        goTerms.append(list(data_go['go-term'][data_go['Gene'].str.match(ii)]))
    
    return goTerms

goTermsInteractors = goTermsFromGenes(uniqueInteractorGenes,data_go)
goTermsQuery = goTermsFromGenes(query,data_go)
    
#%% 
#%% getting and saving 1st layer children of GO terms to file (ONLY RUN ONCE AS IT TAKES FOREVER)
# probably using a dict rather than 2 lists makes more sens

# save all 1st layer children for a gene or save all 1st layer children for a GO terms? latter would take more time later but saves runtime here...
# may also just do this for all 6000 genes so we are done...


def childrenFromGoTerms(genes,goTerms):
    '''
    Returns a dataframe of all 1st layer children per genes based on supplied go terms
    Requires genes and goTerms to be matching & ordered...
    '''
    childrenGoTerms = defaultdict(dict)
    jj = 0    
        
    for ii in genes: #this loop takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?
        tmp = [] #clear tmp variable.
        print('Progress: gene ' + str(jj+1) + '/' + str(len(genes)))
        queryGoTerm = goTerms[jj]
        for kk in queryGoTerm:
            tmp = tmp + getChildrenGoTerm(kk) #This method results in the double checking of a lot of GO terms... Maybe better to get a list of 
            
        childrenGoTerms[ii] = tmp  #does this method overwrite the old values? probably... 
        jj = jj+1
        
    df=pd.DataFrame([childrenGoTerms]).T
    return df

genes = uniqueInteractorGenes[0:2] + uniqueInteractorGenes[-3:-1]

df = childrenFromGoTerms(genes,goTermsInteractors)
df2 = childrenFromGoTerms(query,goTermsQuery)

#%% saving dfs
df.to_excel(r'../data/1stLayerGO_INT_BEM1_BEM2.xlsx', index = True)
df2.to_excel(r'../data/1stLayerGO_BEM1_BEM2.xlsx', index = True)

#%% for each go term of each interactor gene, get its 1st layer children
# for each interactor, find the overlapping set of 1st layer GO terms with our query gene

data_childrenGoTerms=pd.read_excel('../data/1stLayerGO_INT_BEM1_BEM2.xlsx')
# data_childrenGoTerms.rename(columns={"unnamed: 0": "Gene", "0": "GoTerms"})
data_childrenGoTerms.columns = ['Gene','ChildGoTerms']

#%% Fix this so it finds the overlap for query gene children rather than current childrenGoTermsQuery
GOQuery = 'establishment or maintenance of cell polarity'
print('using ' + GOQuery + ' as parent query GO term')
childrenGoTermsQuery = getChildrenGoTerm(GOQuery) #should we do this for every query gene? check plan on github
commonChildGoTerms = defaultdict(dict)
jj = 0

uniqueInteractorGenesTest = uniqueInteractorGenes

for ii in uniqueInteractorGenesTest: #this loop takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?
    print(jj)
    tmp =  list(data_childrenGoTerms['ChildGoTerms'][data_childrenGoTerms['Gene']==ii])
    childrenGoTermsInteractor = tmp[0].split("', '") # bit of a hack... also the first and last entry contain a bracket
    childrenGoTermsInteractor[0].replace("[","")
    childrenGoTermsInteractor[-1].replace("]","")
    tempCommonChildGoTerms = childrenGoTermsQuery + childrenGoTermsInteractor
    commonChildGoTerms[ii] = set([x for x in tempCommonChildGoTerms if tempCommonChildGoTerms.count(x) > 1]) #finds all common go terms - empty so far, does it work properly?
    jj = jj+1
# GOTermsInteractors = getGOTermsGene() 



