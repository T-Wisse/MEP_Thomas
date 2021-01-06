# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 22:00:23 2020

@author: Thomas
"""

from intermine.webservice import Service
import numpy as np
import scipy.io
import seaborn as sns
from scipy import stats, optimize, interpolate
import pandas as pd
from collections import defaultdict 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from scipy import stats
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os, fnmatch

import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
#rel_path_domains="datasets/proteins-domains-from-Pfam.xlsx"

# abs_file_path_domains = os.path.join(script_dir, rel_path_domains)

# os.chdir('../') #<-- for binder os.chdir('../')
# my_path_domains=abs_file_path_domains
#data_domains=pd.read_excel(my_path_domains,header=0,index_col='Unnamed: 0')
#data_domains=data_domains.dropna()

#data_domains.head()

pathways=[ 'galactose degradation','phospholipid biosynthesis', "ethanol degradation"]
biological_processes_slim=["cell budding","lipid binding","cytokinesis"]
biological_processes_goterm=["cell morphogenesis involved in conjugation with cellular fusion","budding cell apical bud growth","establishment of cell polarity","cytoskeleton organization"]
def from_go_to_genes(go,label):
    #label=["GOTerm" or "GOSlimTerm"]
    service = Service('https://yeastmine.yeastgenome.org/yeastmine/service', token = 's1Z7G6d1KdH91d28kd70')
    query = service.new_query("Gene")
    query.add_constraint("goAnnotation.ontologyTerm.parents", label)
    query.add_view(
        "symbol", "goAnnotation.evidence.code.annotType",
        "goAnnotation.ontologyTerm.parents.name"
    )
    query.add_constraint("goAnnotation.qualifier", "!=", "NOT", code = "C")
    query.add_constraint("goAnnotation.qualifier", "IS NULL", code = "D")
    query.add_constraint("goAnnotation.evidence.code.annotType", "=", "manually curated", code = "F")
    query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", go, code = "G")
    query.set_logic("(C or D) and F and G")

    data_toy=defaultdict(dict)

    for row,counter in zip(query.rows(),np.arange(0,len(query.rows()))):

        data_toy['gene-name'][counter]=row["symbol"]
        data_toy['evidence'][counter]=row["goAnnotation.evidence.code.annotType"]
        data_toy['annotation'][counter]=row["goAnnotation.ontologyTerm.parents.name"]


    data_toy_pd=pd.DataFrame(data_toy)
    data=data_toy_pd.drop_duplicates()

    data.index=np.arange(0,len(data))
    return data

#%%

go=biological_processes_goterm[2]
data=from_go_to_genes(go,label='GOTerm')

#%%

plt.plot([1,2,3,4],[1,2,3,4])
plt.ylabel('number of common interactors')
