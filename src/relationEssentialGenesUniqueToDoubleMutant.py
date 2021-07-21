# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:11:22 2021

@author: Thomas
"""


import os
import sys

fileDir = os.path.dirname(os.path.abspath(__file__)) # get path to this file
os.chdir(fileDir) # Change working directory to current file
sys.path.insert(0, fileDir+'\python_modules') #hack to add module to path. otherwise it won't be found. Should find a better way to do this.


from modulesRelationEssentialGenes import getGeneticInteractionsFromGene
from collections import defaultdict
import pandas as pd

#%% get data

essentialGenes_DD = pd.read_excel('../data/essentialGenesDPL1_PSD2_Wessel.xlsx') #essential genes of the double delete. In this case DPL1 & PSD2
essentialGenes_DD = essentialGenes_DD.drop('note','columns')

deletedGenes = ['DPL1','PSD2']

#%% get genes that have a synthetic lethal interaction with query gene
queryGene = 'bem1'

geneticInteractions = getGeneticInteractionsFromGene(queryGene) #NOTE: Returns some 'none' values as interacting genes. WHY?
essentialGenes = geneticInteractions['interactingGene'][geneticInteractions['interactionType'] == 'Synthetic Lethality']

