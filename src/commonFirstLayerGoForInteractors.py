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
import seaborn as sns
from datetime import date

fileDir = os.path.dirname(os.path.abspath(__file__)) # get path to this file
os.chdir(fileDir) # Change working directory to current file

import sys
sys.path.insert(0, fileDir+'\python_modules') #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

from modulesT import common_interactors_T
from modulesT import childrenFromGoTerms
from modulesT import goTermsFromGenes

# from python_modules.module_common_measures import common_partners
#%% User settings
save = False #Set if data should be saved to excel
if save == True:
    dateToday = str(date.today())

#%% getting the data

# data_go = pd.read_excel('../data/slim-goterms-filtered-data.xlsx') #Switched to get data from Yeastmine __> currently unused # Get pre-downloaded go slim data for many genes
# data_go = data_go.dropna()

# data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type'] #Add descriptive headers to columns
# data_go = data_go.drop(['gene-id', 'go-aspect','go-id','feature-type'], axis = 1) #drop unused columns

data=pd.read_excel('../data/data-BioGrid-Yeast.xlsx') #Get pre-downloaded genetic interaction data
data = data.drop(['paper-source'], axis = 1) #drop unused column. SHOULD ALSO DROP INTERACTION THAT ARE NOT SL, NG or PG! Currently this is done before plotting.
#%% query

query = ['BEM1', 'BEM3'] #testing query

# query = pd.read_excel('../data/genesCellPolarity_SGD_amigo.xlsx') #gene set to investigate. Example: cell polarity genes
# query.columns=['gene']
# query = list(set(query.gene)) #get unique entries (remove duplicates in query list)


#%% (TESTING) Reduce dataset to only query genes

# data_go_temp = data_go[data_go['Gene'].isin(query)] #take only GO data of genes which are in the query
# data_temp = data[data['gene-query-name'].isin(query)] #take only interaction data of genes which are in the query. Currently unused


#%% Calling the function common_interactors
# start = timer()
commonInteractorData=common_interactors_T(query,data) #find all common interactors of each of the query genes (with all yeast genes). Based on supplied interaction data
# end = timer()
# print('Duration of determining interactors of query genes function = ' + str(round(end-start)) + 's') #print duration of function common_interactors_T
uniqueInteractorGenes = list(set(commonInteractorData['interactorGene'])) #Get unique interactors (remove duplicates to reduce the list of interactor genes to unique genes)



#%% get go term data for genes (Set up for testing: only take the first x genes)

goTermsInteractors = goTermsFromGenes(uniqueInteractorGenes[0:4])
goTermsQuery = goTermsFromGenes(query)
# goTermsInteractors.columns=['go-term']
# goTermsQuery.columns=['go-term']


#%% getting 1st layer children of GO terms to file (Set up for testing: only take the first x genes)

# save all 1st layer children for a gene or save all 1st layer children for a GO terms? latter would take more time later but saves runtime here...
# may also just do this for all 6000 genes so we are done...
# this fucntion takes essentially forever --> find a better way. maybe download all children for GO terms so we don't have to use yeastmine everytime?

fLCInteractors = childrenFromGoTerms(uniqueInteractorGenes[0:4],goTermsInteractors) #fLC = first layer children. 
fLCQuery = childrenFromGoTerms(query,goTermsQuery)

#%% return unique interactor list to per interacting genes

commonInteractorSetData = defaultdict(dict)

for ii in commonInteractorData['commonInteractors'].to_dict():
    interactorGene = ii.split("-",1)[1]
    commonInteractorSetData[ii] = fLCInteractors[interactorGene]


#%% saving 1st layer children data to excel (using dataframes, if this turns out to be an issue switch to csv package maybe)

if save == True:
    df=pd.DataFrame([fLCInteractors]).T
    df2 = pd.DataFrame([fLCQuery]).T
    
    df.to_excel(r'../data/' + dateToday + '_1stLayerGO_INT_BEM1_BEM3_testing.xlsx', index = True)
    df2.to_excel(r'../data/' + dateToday + '_1stLayerGO_BEM1_BEM3_testing.xlsx', index = True)

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

for goQuery in fLCQuery:
    
    #loop over list of all of all interactors for gene goQuery
    interactingGenes = list(data['gene-target-name'][data['gene-query-name']==goQuery])
    
    for ii in interactingGenes:
        
        tempCommonChildGoTerms = list(set(fLCQuery[goQuery])) + list(set(fLCInteractors[ii])) #Combine unique entries of first layer children into a list
        commonChildGoTerms[goQuery +'-'+ ii] = set([x for x in tempCommonChildGoTerms if tempCommonChildGoTerms.count(x) > 1]) #finds all common go terms - empty so far, does it work properly?
        commonChildGoTermsFraction[goQuery +'-'+ ii] = len(commonChildGoTerms[goQuery +'-'+ ii])/len(fLCQuery[goQuery]) #calculate fraction of common 1st layer children       
        
    jj = jj+1  
    print('Progress: '+ str(jj) + '/'+ str(len(fLCQuery)))

        
#%% saving 1st layer child data to excel

if save == True:
    dfOut=pd.DataFrame([commonChildGoTerms]).T
    dfOut.to_excel(r'../data/' + dateToday+ '_common1stLayerGO_BEM1_BEM3_testing.xlsx', index = True)

#%% reading data saved above
# testing=pd.read_excel('../data/common1stLayerGO_BEM1_BEM3_testing.xlsx')
# testing.columns = ['Genes','ChildGoTerms']


#%% consolidate data for plotting

# interactionType = dict(zip(data['gene-query-name']+'-'+ data['gene-target-name'], data['interaction-type'])) #transform dataframe to dict with keys synthesised from two colums

dataCommonGo = commonInteractorData
dataCommonGo['commonChildGoTermsFraction']= dataCommonGo.index.to_series().map(commonChildGoTermsFraction) #map by index to ensure correct order
dataCommonGo['commonInteractorCount'] = list(map(len,list(commonInteractorData['commonInteractors']))) #superfluous or not?
#dataCommonGo['interactionType'] = commonInteractorData['interactionType']

#%% select which data to plot
gene = 'BEM1'
plotData = dataCommonGo[dataCommonGo['queryGene']==gene]


plotDataReduced = plotData[np.array(plotData['interactionType'] == 'Negative Genetic') | \
                           np.array(plotData['interactionType'] == 'Positive Genetic') | \
                           np.array(plotData['interactionType'] == 'Synthetic Lethality')] #only select specific interaction types

# keys = list(commonChildGoTerms.keys())
# queryKeys = [i for i in keys if gene in i] #Get a list all keys containing the string <gene> 
# fractionCommonGO = list(map(commonChildGoTermsFraction.get, queryKeys)) #get all fractions for query keys


#%% Plotting
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(plotDataReduced,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType', \
                   hue_order=['Negative Genetic','Positive Genetic','Synthetic Lethality'], \
                   diag_kind="hist", diag_kws = {'bins':int(np.ceil(np.sqrt(plotDataReduced.shape[0])))},corner=True)
# plot=sns.pairplot(plotDataReduced,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType',
                  # corner=True, diag_kws = {'bw' : 10, 'kernel' : 'tri'})
plot.fig.suptitle(gene)


#%% Plotting results: Combined plot

# gene = 'BEM1'
# bins = int(np.ceil(np.sqrt(tmp.shape[0])))
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(plotDataReduced,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType', \
                    hue_order=['Negative Genetic','Positive Genetic','Synthetic Lethality'], \
                    diag_kind="hist",diag_kws = {'bins':int(np.ceil(np.sqrt(plotDataReduced.shape[0])))},corner=True)
# plot=sns.pairplot(plotDataReduced,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType', \
#                    hue_order=['Negative Genetic','Positive Genetic','Synthetic Lethality'])
# plt.title(query[0])
plot.fig.suptitle('BEM1',y=1.08)
# # plot.savefig('../output_images/Testing/common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)
