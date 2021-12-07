# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:57:47 2020

Working version to find all 1st layer children GO terms for interacting genes and compare this with common interactor count 

@author: Thomas
"""



import numpy as np
from collections import defaultdict 
import pandas as pd
import os
import seaborn as sns
import scipy as sc
import matplotlib.pyplot as plt
from datetime import date


fileDir = os.path.dirname(os.path.abspath(__file__)) # get path to this file
os.chdir(fileDir) # Change working directory to current file

import sys
sys.path.insert(0, fileDir+'\python_modules') #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

from modulesT import common_interactors_T
from modulesT import childrenFromGoTerms
from modulesT import goTermsFromGenes

#%% User settings
save = False #Set if data should be saved to excel
savePlots = True #Set if plots should be saved

if save|savePlots == True:
    dateToday = str(date.today())

#%% getting locally saved data

# data_go = pd.read_excel('../data/slim-goterms-filtered-data.xlsx') #Switched to get data from Yeastmine __> currently unused # Get pre-downloaded go slim data for many genes
# data_go = data_go.dropna()

# data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type'] #Add descriptive headers to columns
# data_go = data_go.drop(['gene-id', 'go-aspect','go-id','feature-type'], axis = 1) #drop unused columns

data=pd.read_excel('../data/data-BioGrid-Yeast.xlsx') #Get pre-downloaded genetic interaction data
data = data.drop(['paper-source'], axis = 1) #drop unused column. NOTE: Currently unused interactions are dropped before plotting. Doing it here instead may shorten run times
#%% Set query gene(s)

# query = ['BEM1', 'BEM3'] #testing query
query = ['DPL1','PSD2']

# query = pd.read_excel('../data/genesCellPolarity_SGD_amigo.xlsx') #gene set to investigate. Example: cell polarity genes
# query.columns=['gene']
# query = list(set(query.gene)) #get unique entries (remove duplicates in query list)



#%% Calling the function common_interactors

# Possible investigation idea: reduce the interaction gene set to just those contained in the query

commonInteractorData=common_interactors_T(query,data) #find all common interactors of each of the query genes (with all yeast genes). Based on supplied interaction data
uniqueInteractorGenes = list(set(commonInteractorData['interactorGene'])) #Get unique interactors (remove duplicates to reduce the list of interactor genes to unique genes)


#%% get go term data for genes (Set up for testing: only take the first x genes)

goTermsInteractors = goTermsFromGenes(uniqueInteractorGenes)
goTermsQuery = goTermsFromGenes(query)


#%% getting 1st layer children of GO terms to file (Set up for testing: only take the first x genes)

# this fucntion takes essentially forever. Maybe download all children for GO terms so we don't have to use yeastmine everytime?
# save all 1st layer children for a gene or save all 1st layer children for all GO terms? latter would probably take too much time later but once done would save a lot of runtime...
# may also just do this for all 6000 genes so we are done


fLCInteractors = childrenFromGoTerms(uniqueInteractorGenes,goTermsInteractors) #fLC = first layer children. 
fLCQuery = childrenFromGoTerms(query,goTermsQuery)

#%% return unique interactor list to per interacting genes

commonInteractorSetData = defaultdict(dict)

for ii in commonInteractorData['commonInteractors'].to_dict():
    commonInteractorSetData[ii] = fLCInteractors[ii.split("-",1)[1]]


#%% saving 1st layer children data to excel (using dataframes, if this turns out to be an issue maybe we can switch to csv package)
#np.savetxt(r'c:\data\np.txt', df.values, fmt='%d') ALTERNATIVE: SAVE TO TXT

if save == True:
    df=pd.DataFrame([fLCInteractors]).T
    df2 = pd.DataFrame([fLCQuery]).T
    
    df.to_excel(r'../data/' + dateToday + '_1stLayerGO_INT_' + query[0] + '_testing.xlsx', index = True)
    df2.to_excel(r'../data/' + dateToday + '_1stLayerGO_' + query[0] + '_testing.xlsx', index = True)

#%% load data saved above

# data_childrenGoTermsInt=pd.read_excel('../data/1stLayerGO_INT_BEM1_BEM3.xlsx')
# data_childrenGoTermsInt.columns = ['Gene','ChildGoTerms']

# data_childrenGoTermsQuery=pd.read_excel('../data/1stLayerGO_BEM1_BEM3.xlsx')
# data_childrenGoTermsQuery.columns = ['Gene','ChildGoTerms']


#%% Find the overlap of all 1st layer children of the GO Terms corresponding to the query gene(s) and its/their interactors

commonChildGoTerms = defaultdict(dict)
commonChildGoTermsFraction = defaultdict(dict)

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
    dfOut=pd.DataFrame([commonChildGoTerms]).T #Something is wrong with this. Found a set() entry in the excel. Maybe also making list of list?
    dfOut.to_excel(r'../data/' + dateToday+ '_common1stLayerGO_' + query[0] + '_testing.xlsx', index = True) 

#%% reading data saved above
# commonChildGoTerms_loaded=pd.read_excel('../data/2021-05-27_common1stLayerGO_' + query[0] + '_testing.xlsx')
# commonChildGoTerms_loaded.columns = ['Genes','ChildGoTerms']


#%% consolidate data for plotting
# interactionType = dict(zip(data['gene-query-name']+'-'+ data['gene-target-name'], data['interaction-type'])) #transform dataframe to dict with keys synthesised from two colums

# commonInteractorData = commonInteractorData.dropna()

dataCommonGo = commonInteractorData
dataCommonGo['commonChildGoTermsFraction']= dataCommonGo.index.to_series().map(commonChildGoTermsFraction) #map by index to ensure correct order

for ii in list(commonInteractorData.index.values):
    print(ii)
    dataCommonGo.loc[ii,'commonInteractorCount'] = len(commonInteractorData.loc[ii,'commonInteractors'])
    
#%% Saving final data
if save == True:
    dataCommonGo.to_excel(r'../data/' + dateToday+ '_common1stLayerGO_BEM1_BEM3_dataCommonGo.xlsx', index = True) 


#%% Plotting 

def getPlotData(gene,dataCommonGo): #interaction data must have a 'gene query name' and a 'gene target name' column
    '''
   
    '''
    # d2 = defaultdict(dict)
   
    plotData = dataCommonGo[dataCommonGo['queryGene']==gene]
    
    plotData = plotData[np.array(plotData['interactionType'] == 'Negative Genetic') | \
                               np.array(plotData['interactionType'] == 'Positive Genetic') | \
                               np.array(plotData['interactionType'] == 'Synthetic Lethality')] #only select specific interaction types
    return plotData


# gene = 'BEM1' #Plot data of this gene
gene = query[0]

plotData = getPlotData(gene,dataCommonGo)

# bins = int(np.ceil(np.sqrt(plotData.shape[0])))
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(plotData,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType', \
                    hue_order=['Negative Genetic','Positive Genetic','Synthetic Lethality'], \
                    diag_kind="hist",diag_kws = {'bins':int(np.ceil(np.sqrt(plotData.shape[0])))},corner=True)
# plt.title(gene)
plot.fig.suptitle(gene,y=1.08)
if savePlots == True:
    plot.savefig('../data/images/1stLayerCommonGO_byType_' + gene + '.png',dpi=300,format='png',transparent=True)



#Alternative plotting
# sns.set(style="ticks", color_codes=True)
# plot=sns.pairplot(plotData,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType', \
#                     hue_order=['Negative Genetic','Positive Genetic','Synthetic Lethality'], \
#                     diag_kind="hist", diag_kws = {'bins':int(np.ceil(np.sqrt(plotData.shape[0])))},corner=True)
# plot=sns.pairplot(plotData,vars=['commonChildGoTermsFraction','commonInteractorCount'],hue='interactionType',
#                    corner=True, diag_kws = {'bw' : 10, 'kernel' : 'tri'})
# plot.fig.suptitle(gene)

#%% Comparing means of the three distributions (SL, NG, PG) Determine p value of equal mean...
# must work for different sample counts
# assuming normal distribution...
# assuming independent samples
# trim outliers?
# pValue < 0.05 means the means of the two distributions are not equal.
[t1,pValue_NG_SL] = sc.stats.ttest_ind(plotData['commonInteractorCount'][plotData['interactionType'] == 'Negative Genetic'], plotData['commonInteractorCount'][plotData['interactionType'] == 'Synthetic Lethality'] , equal_var=False)
[t2,pValue_NG_PG] = sc.stats.ttest_ind(plotData['commonInteractorCount'][plotData['interactionType'] == 'Negative Genetic'], plotData['commonInteractorCount'][plotData['interactionType'] == 'Positive Genetic'] , equal_var=False)

#%% Comparing same genetic interaction between different genes. 
gene_2 = query[1] #Plot data of this gene

plotData_2 = getPlotData(gene_2,dataCommonGo)

#CommonChildGoTermsFraction
[t3,pValue_genes_NG] = sc.stats.ttest_ind(plotData['commonChildGoTermsFraction'][plotData['interactionType'] == 'Negative Genetic'], plotData_2['commonChildGoTermsFraction'][plotData_2['interactionType'] == 'Negative Genetic'] , equal_var=False)

