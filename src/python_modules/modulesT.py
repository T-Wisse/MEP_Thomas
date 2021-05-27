# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:58:56 2020
@author: Thomas

Contains functions

goTermsFromGenes    Currently a bit convoluted. May have broken other scripts by removing getGoTermsFromGene! --> update to new function
getChildrenGoTerm   Currently a bit convoluted. May have broken other scripts by removing getGoTermsFromGene! --> update to new function

getPathwayForGeneFunc 
common_interactors_T 
common_go_children 



"""
from __future__ import print_function
from intermine.webservice import Service
from collections import defaultdict 
import numpy as np
import pandas as pd

# from collections import defaultdict 

def goTermsFromGenes(genes): #get go terms for genes from yeastmine
    '''
    Returns a dict of all go terms for supplied genes
    
    '''
    
    def getGoTermsFromGene(gene): #
        '''
        Returns a list of all the GO terms of the supplied gene (name) from yeastmine.
        
        paremeters:
        -------
        gene: string
        '''
      
        from intermine.webservice import Service
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query("Gene")
        query.add_view("symbol", "goAnnotation.ontologyTerm.name")
        query.add_sort_order("Gene.primaryIdentifier", "ASC")
        query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="F")
        query.add_constraint("status", "IS NULL", code="C")
        query.add_constraint("status", "=", "Active", code="B")
        query.add_constraint("Gene", "LOOKUP", gene, code="A")
        query.set_logic("(B or C) and F and A")
        
        GOterms = []
        
        for row in query.rows():
            GOterms.append(row["goAnnotation.ontologyTerm.name"])
        
        GOterms = list(set(GOterms)) #get unique values
    
        return GOterms

    goTerms = defaultdict(dict)
    jj = 0
        
    for ii in genes: #For each gene in the supplied set of genes. #this loop takes essentially forever --> Maybe download all children for GO terms so we don't have to use yeastmine everytime?
        # tmp = [] #clear tmp variable.
        print('Progress: gene ' + str(jj+1) + '/' + str(len(genes)))
        jj = jj+1
        try: 
            goTerms[ii] = getGoTermsFromGene(ii) #This method results in the double checking of a lot of GO terms... Maybe better to get a list of 
        except Exception: #Specific error I want to catch is WebserviceError which is part of the intermine.webservice package but I don't know how
            goTerms[ii] = 'Error'
    # df=pd.DataFrame([goTerms]).T #transform dict to dataframe.
    return goTerms

def childrenFromGoTerms(genes,goTerms): 
    '''
    Returns a dataframe of all 1st layer children per genes based on supplied go terms
   
    '''
    
    def getChildrenGoTerm(GOTerm):
       
        '''
        Returns a list of all the children of the supplied GO term (name) from yeastmine.
        If list is empty, check if a valid GO Term name is supplied
        
        paremeters:
        -------
        GOTerm: str (name)
        '''
    
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query("Gene")
        query.add_constraint("goAnnotation.ontologyTerm.parents", "GOTerm")
        query.add_constraint("goAnnotation.ontologyTerm", "GOTerm")
        query.add_view("goAnnotation.ontologyTerm.name", "goAnnotation.ontologyTerm.identifier")
        query.add_constraint("organism.name", "=", "Saccharomyces cerevisiae", code="B")
        query.add_constraint("status", "IS NULL", code="F")
        query.add_constraint("status", "=", "Active", code="E")
        query.add_constraint("goAnnotation.qualifier", "IS NULL", code="J")
        query.add_constraint("goAnnotation.qualifier", "!=", "NOT", code="I")
        query.add_constraint("goAnnotation.evidence.code.annotType", "=", "manually curated", code="D")
        query.add_constraint("goAnnotation.evidence.code.annotType", "=", "high-throughput", code="C")
        query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", GOTerm, code="A")
        query.set_logic("A and B and (C or D) and (E or F) and (I or J)")
        
        GOChildren = []
        
        for row in query.rows():
            GOChildren.append(row["goAnnotation.ontologyTerm.name"])
            
        GOChildren = list(set(GOChildren)) #get unique values
        
        return GOChildren
    
    childrenGoTerms = defaultdict(dict)
    jj = 0
     
    for ii in genes: #For each gene in the supplied set of genes. 
        tmp = [] #clear tmp variable.
        print('Progress: gene ' + str(jj+1) + '/' + str(len(genes)))
        jj = jj+1
        
        queryGoTerm = goTerms[ii] #Select corresponding set of go terms for the gene
        for kk in queryGoTerm: #For each go term of the gene, find all 1st layer children 
            for attempt in range(2): #Allow for multiple attempts
                try: 
                    tmp = tmp + getChildrenGoTerm(kk) #This method results in the double checking of a lot of GO terms... Maybe better to get a list of 
                except Exception: #catch webservice errors
                    tmp = tmp + ['Error_attempt_'+str(attempt+1)] #first draft of error messages. Probably not very useful like this but eh
                else: #if no exception is throw, break the for loop
                    break
                
        if queryGoTerm in tmp: #check if list contains original (parent) go term and remove it
            tmp.remove(queryGoTerm)
        childrenGoTerms[ii] = tmp
                  
    return childrenGoTerms




def getPathwayForGeneFunc(gene):
        
    
    dictOut = defaultdict()
    
    for ii in gene:
        
        tmp = []        
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query("Gene")
        query.add_view("symbol", "pathways.identifier", "pathways.name")
        query.add_sort_order("Gene.primaryIdentifier", "ASC")
        query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="B")
        query.add_constraint("Gene", "LOOKUP", ii, code="A")
        
        for row in query.rows():
            # dictOut[ii] = (row["symbol"], row["pathways.identifier"], row["pathways.name"])
            tmp.append(row["pathways.identifier"])
        
        dictOut[ii] = tmp
    return dictOut

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return(a_set.intersection(b_set))  
    return()

def common_interactors_T(query,interactionData): #interaction data must have a 'gene query name' and a 'gene target name' column
    '''
    Finds common interactors for all genes in the query based on supplied interaction data
    
    --------
    paremeters:
    query: list of genes
    interactionData: dataframe of genetic interactions, contains gene query, gene target and type of interaction
    --------
    output: dataframe with interacting genes, common interactors and interaction type        
    '''
    # d2 = defaultdict(dict)
    commonInteractors = defaultdict(dict)
    
    
    for queryGene in query: #For each query gene
    
        query1 = interactionData[interactionData['gene-query-name']==queryGene] #Reduce list to only query gene
        interactorList = query1['gene-target-name'].unique() #List genes interacting with query1 gene 
        # interactionTypes = interactionData.loc[interactionData['gene-query-name'] == queryGene, 'interaction-type']
        
        # print(interactorList)
        # print(interactionTypes)

        # print(interactionTypeList)
        for interactorGene in interactorList: #For each interactor 

            query2 = interactionData[interactionData['gene-query-name']==interactorGene]
            interactors = query2['gene-target-name'].unique() #List genes interacting with query2 gene
            tempCommonInt = list(common_member(interactorList, interactors))
                        
            tmpInteractionType = query1[query1['gene-target-name']==interactorGene]['interaction-type'].unique().tolist()
                      
            
            commonInteractors[queryGene+'-'+interactorGene]['interactionType'] = tmpInteractionType[0]
            commonInteractors[queryGene+'-'+interactorGene]['queryGene'] = queryGene
            commonInteractors[queryGene+'-'+interactorGene]['interactorGene'] = interactorGene
            commonInteractors[queryGene+'-'+interactorGene]['commonInteractors'] = tempCommonInt
            # tmpInteraction = query1[]
            
            # commonInteractors[queryGene+'-'+interactorGene]['interactionType'] = query1['gene-target-name'] ==
        
            # for interactorName  in interactors:
            #     if interactorName in interactorList: # if a gene interactor of the query1 is in interactors of query 2 
            #         commonInteractors.append(interactorName) #add gene to list
      
    df=pd.DataFrame(commonInteractors).T #interactors in df has to be hashable list (no list of lists)
    
    return df


def common_go_children(data_go,data_common_partners):
    """"
    function that computes the common go terms or interactors genes
    input: data_go= dataframe with all go terms per gene
    data_common_partners= dataframe output of the function common_partners
    
    """
    d3=defaultdict(dict)
    query=np.unique(np.array(data_common_partners['query'])) # Finds all query genes listed with a common partner
    # big for loop for each gene analyzed in common partners
    for i in np.arange(0,len(query)):
        partners=data_common_partners[data_common_partners['query']==query[i]]['names of genes']

        d2=defaultdict(dict)
        
        # print(partners)
        # print(query)
        for genes in partners:
            d2[genes]['query']=query[i]
            d2[genes]['names of genes']=genes

            tmp=data_go[data_go['Gene']==query[i]]['go-term'].tolist()
            tmp=np.unique(tmp).tolist()

            

            tmp2=data_go[data_go['Gene']==genes]['go-term'].tolist()
            tmp2=np.unique(tmp2).tolist()

                        

       
            d2[genes]['common-go-terms']=np.intersect1d(tmp,tmp2)
            if len(tmp)==0:
                d2[genes]['fraction-of-common-go']=0
            else:
                d2[genes]['fraction-of-common-go']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100
                
      
        d3.update(d2)
        
    df=pd.DataFrame(d3).T
    df_sorted=df.sort_values(by='fraction-of-common-go',ascending=False)

    return df_sorted
