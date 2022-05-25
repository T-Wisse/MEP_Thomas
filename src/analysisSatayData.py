# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:00:41 2021

@author: Thomas
"""

import pandas as pd
from collections import defaultdict 

#%%
######## vOLCANO PLOTS ######## ???

dataDir = 'D:/SATAY/SATAY_Data/SATAY_Analysis/' # Directory containing processed SATAY data

strains = ['yTW001_4','yTW001_6','yLIC137_7','yLIC137_8']
# strains = ['yLIC137_7'] #Run a single strain for testing purposes

dataStrains = pd.DataFrame() #initialize dataframe

for ii in strains:
    datafile = dataDir + ii + '/' + ii + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt' # Select data file
    
    rawData = pd.read_csv(datafile, sep="\t", index_col=0) #read txt file, use tab as seperator, column 0 (gene name) becomes the index
    rawData.columns = rawData.columns.str.title() # Capitalize each word
    rawData.columns = rawData.columns.str.replace(' ', '') #Remove whitespaces in column headers
    rawData.drop("ADE2",axis=0,inplace=True) #remove the ADE2 entry
    
    minReadCount = 10 #set a minumum number of reads
    data = rawData[rawData.NumberOfReadsPerGene >= minReadCount] #Drop entries with read count under minimum (?)
    
    data['readsPerTransposon'] = data['NumberOfReadsPerGene']/data['NumberOfTransposonsPerGene']
    # transposonsPerGene = data['readsPerTransposon'].to_dict() #Get reads/transposon as dict for calculations
    
    
    dataStrains[ii] = data['readsPerTransposon'] #combine reads/transposon data for each strain
    
    
#%% differences between biological (?) replicates 




#%% differences between technical replicates


# df.to_dict()


#%% Call transposonmapper functions
# OUTDATED!
# import os, sys
# import matplotlib.pyplot as plt
# # from D:\Users\Thomas\Studie\MEP\Transposonmapper\transposonmapper\statistics.py import volcano
# # from transposonmapper.statistics.volcano_helpers import apply_stats,  info_from_datasets, make_datafile
# from transposonmapper import volcano

# a = 'yTW001_4'
# b = 'yTW001_6'

# path_a =  datafile = dataDir + a# Select data file a
# filelist_a = a + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt' 
# path_b = datafile = dataDir + b # Select data file b
# filelist_b = b + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt'

# volcano(path_a, filelist_a, path_b, filelist_b, variable='read_per_gene', significance_threshold=0.01, normalize=True, trackgene_list=[], figure_title="test")

#%% Call transposonmapper functions

import os, sys
import matplotlib.pyplot as plt
# from D:\Users\Thomas\Studie\MEP\Transposonmapper\transposonmapper\statistics.py import volcano
# from transposonmapper.statistics.volcano_helpers import apply_stats,  info_from_datasets, make_datafile
from transposonmapper import statistics


a1 = 'yTW001_4'
a2 = 'yTW001_6'
b1 = 'yLIC137_7'
b2 = 'yLIC137_8'

dataDirYTW = dataDir + 'yTW001/'
dataDirYLIC = dataDir + 'yLIC137/'

path_a = dataDirYTW # Select data file a
filelist_a = [dataDirYTW + a1 + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt', dataDirYTW + a2 + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt']
path_b = dataDirYLIC # Select data file b
filelist_b = [dataDirYLIC + b1 + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt', dataDirYLIC + b2 + '_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene.txt']

bonferroniThreshold= 0.05/6600 #bonferoni correction
orginalThreshold = 0.01 #orginally used

threshold = bonferroniThreshold
#WHAT ABOUT FWER? SEEMS THE MOST REASONALBE TO ME. We want to reduce number of errors but don't really care too much if we have false positives or false negatives

dfTnPerGene = statistics.volcano(path_a, filelist_a, path_b, filelist_b, variable='tn_per_gene', significance_threshold=threshold, normalize=True, trackgene_list=[], figure_title="yTW001 vs yLIC137")
dfReadsPerGene = statistics.volcano(path_a, filelist_a, path_b, filelist_b, variable='read_per_gene', significance_threshold=threshold, normalize=True, trackgene_list=[], figure_title="yTW001 vs yLIC137")
dfReadsPerInsrt = statistics.volcano(path_a, filelist_a, path_b, filelist_b, variable='Nreadsperinsrt', significance_threshold=threshold, normalize=True, trackgene_list=[], figure_title="yTW001 vs yLIC137")

#%% apply FDR to determine significance
import numpy as np

alpha = 0.05
m = 6600
k = np.array(range(1,6601))

dfTnPerGene = dfTnPerGene.sort_values('p_value',ascending=False) # Note: p values are retured as -1*np.log10(pval). So we need to do 10^-pvalue to obtain the actual pvalue. Alternatively, we can apply teh same math to the fdr value for comparison
dfReadsPerGene = dfReadsPerGene.sort_values('p_value',ascending=False) # Note: p values are retured as -1*np.log10(pval). So we need to do 10^-pvalue to obtain the actual pvalue. Alternatively, we can apply teh same math to the fdr value for comparison
dfReadsPerInsrt = dfReadsPerInsrt.sort_values('p_value',ascending=False) # Note: p values are retured as -1*np.log10(pval). So we need to do 10^-pvalue to obtain the actual pvalue. Alternatively, we can apply teh same math to the fdr value for comparison

sTnPerGene = np.array(dfTnPerGene['p_value'])
sReadsPerGene = np.array(dfReadsPerGene['p_value'])
sReadsPerInsrt = np.array(dfReadsPerInsrt['p_value'])

FDRvalues = k * alpha / m

pValuesTnPerGene = np.power(10,-sTnPerGene)
pValuesReadsPerGene = np.power(10,-sReadsPerGene)
pValuesReadsPerInsrt = np.power(10,-sReadsPerInsrt)

significanceFDRTnPerGene = (pValuesTnPerGene < FDRvalues)
significanceFDRReadsPerGene = (pValuesReadsPerGene < FDRvalues)
significanceFDRReadsPerInsrt = (pValuesReadsPerInsrt < FDRvalues)


#%%
## Importing the required python libraries 
import os, sys
import warnings
import timeit
import numpy as np
import pandas as pd 
import pkg_resources
from transposonmapper.processing.clean_bedwigfiles import cleanfiles


#%%
######## GENOME PROFILE PLOTS ########

######## Lets save the wig and bed files as variables to clean them and call the function#####################

wig_files=[]
bed_files=[]

sample = "yLIC137_8"
data_dir="C:/Users/Thomas/Desktop/Uni/SATAY_Data/rawDataFiles/" #Folder contains all the un-edited output files of the transposonmapper for every sample

wig_files = data_dir + sample + '_merged_cleaned_forward_reads_trimmed.sorted.bam.wig'
bed_files = data_dir + sample + '_merged_cleaned_forward_reads_trimmed.sorted.bam.bed'



#%% Visualize the insertions and reads per gene throughout the genome
 ## Import the function
 
from transposonmapper.processing.transposonread_profileplot_genome import profile_genome

#### vizualization #####
bed_file= bed_files # example for the 1st file 
bar_width=None
savefig=False

# Variable may be "reads" or "transposons" #NOTE: READS GIVES A HUGE PEAK AT 1 POINT. SOMETHING IS EITHER WRONG OR WE NEED TO CUT THIS PEAK FROM THE FIGURE FOR BETTER SCALING
# Outliers for gene #233
# Edit D:\Anaconda\envs\satay\lib\site-packages\transposonmapper\plotting\profile_genome_plot.py for editing the plot
# !!! HAVE TO EDIT THE PLOT TITLE BY HAND IN THAT FILE! !!!
profileTransposons=profile_genome(bed_file=bed_file, variable="transposons", bar_width=bar_width, savefig=savefig,showfig=True)
profileReads=profile_genome(bed_file=bed_file, variable="reads", bar_width=bar_width, savefig=savefig,showfig=True)
