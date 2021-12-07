# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 09:28:15 2021

@author: Thomas & Leila (python-modules-for-bioinformatic-analyses)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gseapy as gp 
from gseapy.plot import barplot, dotplot

#%% select data

data_scores=pd.read_excel('../data/scores_wt_nrp1_satay.xlsx',index_col='Unnamed: 0')

#%% Constructing the network matrix to cytoscape
columns=['source_node','target_node','score','type']

data_network=pd.DataFrame(columns=columns)

data_network['source_node'] =  np.ones(shape=(len(data_scores)))
ref=data_scores.loc['HO','score']
data_scores.loc[:,'score_HO']=data_scores.loc[:, 'score']/ ref
col='score'
for i in np.arange(0,len(data_network)):
    data_network.loc[i,'source_node']= 'nrp1'
    j=data_scores.index[i]
    if data_scores.loc[j,col]> data_scores[col].mean()+3.5*data_scores[col].std():
        data_network.loc[i,'target_node']=j
        data_network.loc[i,'score']=data_scores.loc[j,col]
        data_network.loc[i,'type']='PG'
    elif data_scores.loc[j,col]< data_scores[col].mean()-3.5*data_scores[col].std():
        data_network.loc[i,'target_node']=j
        data_network.loc[i,'score']=data_scores.loc[j,col]
        data_network.loc[i,'type']='SL'
    else:
        
        data_network.loc[i,'target_node']=j
        data_network.loc[i,'score']=data_scores.loc[j,col]
        data_network.loc[i,'type']='Neutral'
#%%
data_sl=data_network[data_network.loc[:,'type']=='SL']
data_pi=data_network[data_network.loc[:,'type']=='PG']

data_neutral=data_network[data_network.loc[:,'type']=='Neutral']

#%%
gene_list=data_sl['target_node'].squeeze().str.strip().tolist()


# gene_list=['CDC24','CDC42', 'NRP1', 'BEM1', 'BEM3', 'BEM2', 'RDI1','CLA4']
gene_list=['ICE2','PSD1', 'DGK1', 'RPL20A', 'GPP1', 'OPI8', 
           'SRB5','VMS1','IBA57','RIC1','END3','YPT6','UMP1',
           'SPB1','SNF7','VPS27','YME1','RPL21A','HIP1', 'TIF2',
           'YNL019C','CLB2']
#%% datasets
yeast = gp.get_library_name(database='Yeast')
sets=[yeast[2],yeast[5],yeast[8] ] #['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018']
#%% enrichment 
i=1
enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=sets[i],
                 organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='datasets/enrich-analysis/TEST',
                 # no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
               )
#%%
results=enr.results

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

plot=ax.scatter(results['Term'], results['Combined Score'], c=results['P-value']
                ,cmap='viridis',s=200,alpha=0.7)
ax.set_ylabel('Combined score')
ax.set_xticklabels(results['Term'], rotation=45, ha='right')
ax.set_title('nrp1_sl_enrichment_satay'+sets[i])
cbar=fig.colorbar(plot, ax=ax)
cbar.ax.set_title("P-value")


#%%
data_pi.to_csv('nrp1-network-pi-genes.csv')
data_sl.to_csv('nrp1-network-sl-genes.csv')
#%%
#a_list = ["abc", "def", "ghi"]
a_list=data_low.values
textfile = open("low_values_genes.txt", "w")
for element in a_list:
    textfile.write(element + "\n")
textfile.close()