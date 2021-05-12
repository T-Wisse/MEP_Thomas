# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:57:47 2020

Working version to find all 1st layer children GO terms for interacting genes

How do we want to save this? check if data types are logical and useful.
Manually check if the output is correct...

@author: Thomas
"""

from timeit import default_timer as timer

# import seaborn as sns
# import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict 
import pandas as pd
import os

fileDir = os.path.dirname(os.path.abspath(__file__)) # get path to this file
os.chdir(fileDir) # Change working directory to current file
modulePath = fileDir+'\python_modules'

import sys
sys.path.insert(0, modulePath) #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

from modulesT import common_interactors_T
from modulesT import common_go_children
from modulesT import getChildrenGoTerm
from modulesT import getGoTermsFromGene

# from python_modules.module_common_measures import common_partners
#%% User settings
save = False #Set if data should be saved to excel

#%% getting the data

data_go = pd.read_excel('../data/slim-goterms-filtered-data.xlsx') #Switch to get data from Yeastmine? Get pre-downloaded go slim data for many genes
data_go = data_go.dropna()

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type'] #Add descriptive headers to columns
data_go = data_go.drop(['gene-id', 'go-aspect','go-id','feature-type'], axis = 1) #drop unused columns

data=pd.read_excel('../data/data-BioGrid-Yeast.xlsx') #Get pre-downloaded genetic interaction data
data = data.drop(['paper-source'], axis = 1) #drop unused column
#%% query

query = ['BEM1', 'BEM3'] #testing query

# query = pd.read_excel('../data/genesCellPolarity_SGD_amigo.xlsx') #gene set to investigate. Example: cell polarity genes
# query.columns=['gene']
# query = list(set(query.gene)) #get unique entries (remove duplicates in query list)


#%% (TESTING) Reduce dataset to only query genes

# data_go_temp = data_go[data_go['Gene'].isin(query)] #take only GO data of genes which are in the query
# data_temp = data[data['gene-query-name'].isin(query)] #take only interaction data of genes which are in the query. Currently unused


#%% Calling the function common_interactors
start = timer()
commonInteractorData=common_interactors_T(query,data) #find all common interactors of each of the query genes (with all yeast genes). Based on supplied interaction data
end = timer()
print('Duration of determining interactors of query genes function = ' + str(round(end-start)) + 's') #print duration of function common_interactors_T

#common_go=common_go(data_go=data_go,data_common_partners=common_partners_data)

#%% get all unqiue interactor genes
uniqueInteractorGenes = list(set(commonInteractorData['interactorGene'])) # remove duplicates to reduce the list of interactor genes to unique genes
# uniqueInteractorGenes = commonInteractorData.interactorGene.unique() #other method?

#%% (testing) get go term data for genes 

def goTermsFromGenes(genes): #get go terms for genes from yeastmine
    '''
    Returns a dict of all go terms for supplied genes
    
    '''
    goTerms = defaultdict(dict)
    jj = 0
        
    for ii in genes: #For each gene in the supplied set of genes. #this loop takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?
        # tmp = [] #clear tmp variable.
        print('Progress: gene ' + str(jj+1) + '/' + str(len(genes)))
        jj = jj+1
        
        goTerms[ii] = getGoTermsFromGene(ii) #This method results in the double checking of a lot of GO terms... Maybe better to get a list of 
        
    # df=pd.DataFrame([goTerms]).T #transform dict to dataframe.
    return goTerms

goTermsInteractors = goTermsFromGenes(uniqueInteractorGenes[0:10])
goTermsQuery = goTermsFromGenes(query)
# goTermsInteractors.columns=['go-term']
# goTermsQuery.columns=['go-term']


#%% getting 1st layer children of GO terms to file 

# save all 1st layer children for a gene or save all 1st layer children for a GO terms? latter would take more time later but saves runtime here...
# may also just do this for all 6000 genes so we are done...


def childrenFromGoTerms(genes,goTerms): 
    '''
    Returns a dataframe of all 1st layer children per genes based on supplied go terms
   
    '''
    childrenGoTerms = defaultdict(dict)
    jj = 0
     
    for ii in genes: #For each gene in the supplied set of genes. #this loop takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?
        tmp = [] #clear tmp variable.
        print('Progress: gene ' + str(jj+1) + '/' + str(len(genes)))
        jj = jj+1
        
        queryGoTerm = goTerms[ii] #Select corresponding set of go terms for the gene
        for kk in queryGoTerm: #For each go term of the gene, find all 1st layer children 
            tmp = tmp + getChildrenGoTerm(kk) #This method results in the double checking of a lot of GO terms... Maybe better to get a list of 
        
        if queryGoTerm in tmp: #check if list contains original (parent) go term and remove it
            tmp.remove(queryGoTerm)
        childrenGoTerms[ii] = tmp
                  
    return childrenGoTerms

fLCInteractors = childrenFromGoTerms(uniqueInteractorGenes[0:10],goTermsInteractors) #fLC = first layer children
fLCQuery = childrenFromGoTerms(query,goTermsQuery)


#%% saving 1st layer children data to excel (using dataframes, if this turns out to be an issue switch to csv package maybe)

if save == True:
    df=pd.DataFrame([fLCInteractors]).T
    df2 = pd.DataFrame([fLCQuery]).T
    
    df.to_excel(r'../data/1stLayerGO_INT_BEM1_BEM3_testing.xlsx', index = True)
    df2.to_excel(r'../data/1stLayerGO_BEM1_BEM3_testing.xlsx', index = True)

#np.savetxt(r'c:\data\np.txt', df.values, fmt='%d') ALTERNATIVE: SAVE TO TXT. MAY BE WORTH LOOKING INTO


#%% load data saved above

# data_childrenGoTermsInt=pd.read_excel('../data/1stLayerGO_INT_BEM1_BEM3.xlsx')
# data_childrenGoTermsInt.columns = ['Gene','ChildGoTerms']

# data_childrenGoTermsQuery=pd.read_excel('../data/1stLayerGO_BEM1_BEM3.xlsx')
# data_childrenGoTermsQuery.columns = ['Gene','ChildGoTerms']


#%% Find the overlap of all 1st layer children of the GO Terms corresponding to the query gene(s) and its interactors

commonChildGoTerms = defaultdict(dict)
commonChildGoTermsFraction = defaultdict(dict)

genes = list(fLCInteractors.keys())

jj = 0

for goQuery in fLCQuery: #THIS LOOP LOOKS AT ALL QUERY GENES WITH ALL GENES, NOT JUST INTERACTING GENES 
    
    for ii in fLCInteractors:
        
        tempCommonChildGoTerms = list(set(fLCQuery[goQuery])) + list(set(fLCInteractors[ii])) #Combine unique entries of first layer children into a list
        commonChildGoTerms[goQuery +'-'+ ii] = set([x for x in tempCommonChildGoTerms if tempCommonChildGoTerms.count(x) > 1]) #finds all common go terms - empty so far, does it work properly?
        commonChildGoTermsFraction[goQuery +'-'+ ii] = len(commonChildGoTerms[goQuery +'-'+ ii])/len(fLCQuery[goQuery]) #calculate fraction of common 1st layer children       
        
    jj = jj+1  
    print('Progress: '+ str(jj) + '/'+ str(len(fLCQuery)))


        
#%%
if save == True:
    dfOut=pd.DataFrame([commonChildGoTerms]).T
    dfOut.to_excel(r'../data/common1stLayerGO_BEM1_BEM3_testing.xlsx', index = True)

#%% reading data saved
# testing=pd.read_excel('../data/common1stLayerGO_test.xlsx')
# testing.columns = ['Genes','ChildGoTerms']

if save == True:
    print('test')
    
#%%
testing = defaultdict(dict)

for goQuery in fLCQuery:    
    for ii in fLCInteractors:
        