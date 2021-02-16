# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:28:26 2021

@author: Thomas
"""

from __future__ import print_function
from intermine.webservice import Service

def getPathwayForGene(gene):
    
    tmp = []
    for ii in gene:

        #!/usr/bin/env python
        
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query(str(ii))
        query.add_view("symbol", "pathways.identifier", "pathways.name")
        query.add_sort_order("Gene.primaryIdentifier", "ASC")
        query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="B")
        query.add_constraint("Gene", "LOOKUP", "fas1", code="A")
        
        for row in query.rows():
            tmp.append(row["symbol"], row["pathways.identifier"], row["pathways.name"])  
        