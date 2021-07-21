# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:15:24 2021

@author: Thomas
"""

from __future__ import print_function
from intermine.webservice import Service
from collections import defaultdict 
import pandas as pd

#%% Get genetic interactions for supplied gene from yeastmine


def getGeneticInteractionsFromGene(gene):
    '''
    gene: string of single gene standard name
    
    '''
    #!/usr/bin/env python
   
    
    service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
    query = service.new_query("Gene")
    query.add_constraint("interactions.participant2", "Gene")
    query.add_view(
        "symbol", "interactions.participant2.symbol",
        "interactions.details.experiment.interactionDetectionMethods.identifier"
    )
    query.add_sort_order("Gene.interactions.details.experiment.interactionDetectionMethods.identifier", "ASC")
    query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="B")
    query.add_constraint("interactions.details.relationshipType", "=", "genetic", code="D")
    query.add_constraint("Gene", "LOOKUP", gene, code="A")
    
    
    outputQueryGene = defaultdict(dict)
    outputInteractingGene = defaultdict(dict)
    outputInteractionType = defaultdict(dict)
    ii = 0
    
    for row in query.rows(): #currently loop only works for a single input gene
        outputQueryGene[ii] = row["symbol"]
        outputInteractingGene[ii] = row["interactions.participant2.symbol"]
        outputInteractionType[ii] = row["interactions.details.experiment.interactionDetectionMethods.identifier"]
        ii = ii+1
        
    outputData = pd.DataFrame([outputQueryGene]).T
    outputData.columns = ['interactingGene']
    outputData['interactingGene']= outputData.index.to_series().map(outputInteractingGene) #map by index to ensure correct order
    outputData['interactionType']= outputData.index.to_series().map(outputInteractionType) #map by index to ensure correct order
    
    return outputData