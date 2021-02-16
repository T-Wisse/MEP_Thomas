# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:58:56 2020
@author: Thomas
"""

import numpy as np
import pandas as pd
from __future__ import print_function
from intermine.webservice import Service
# from collections import defaultdict 

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return(a_set.intersection(b_set))  
    return()

def common_partners_T(query,interactionData): #interaction data must have a 'gene query name' and a 'gene target name' column
    '''
    Finds common interactors for all genes in the query based on supplied interaction data
    query: list of genes
    ineractionData: dataframe of genetic interactions, contains gene query, gene target and type of interaction
    '''
    # d2 = defaultdict(dict)
    commonInteractors = []
    
    for queryGene in query: #For each query gene
    
        query1 = interactionData[interactionData['gene-query-name']==queryGene] #Reduce list to only query gene
        interactorList=query1['gene-target-name'].unique() #List genes interacting with query1 gene 
        
        for interactorGene in interactorList:
            
            query2 = interactionData[interactionData['gene-query-name']==interactorGene]
            interactors=query2['gene-target-name'].unique() #List genes interacting with query2 gene
                
            commonInteractors.append(common_member(interactorList, interactors))
            
        
            # for interactorName  in interactors:
            #     if interactorName in interactorList: # if a gene interactor of the query1 is in interactors of query 2 
            #         commonInteractors.append(interactorName) #add gene to list
                    
    return commonInteractors #note, so far will only work properly for 1 gene in query...


# def common_go(data_go,data_common_partners):



def getPathwayForGene(gene):
    
    '''
    Gets all annotated pathways of the query gene(s) from yeastmine
    ------
    parameters:
        gene: a list of one or more query genes. Type: string
    '''
    
    tmp = []
    for ii in gene:

        #!/usr/bin/env python
        
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query(ii)
        query.add_view("symbol", "pathways.identifier", "pathways.name")
        query.add_sort_order("Gene.primaryIdentifier", "ASC")
        query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="B")
        query.add_constraint("Gene", "LOOKUP", "fas1", code="A")
        
        for row in query.rows():
            tmp.append(row["symbol"], row["pathways.identifier"], row["pathways.name"])  
    return tmp