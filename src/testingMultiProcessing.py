# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:33:40 2021

@author: Thomas
"""

path = 'D:/Users/Thomas/Studie/MEP/MEP_Thomas/src/python_modules'
import sys
sys.path.insert(0, path) #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

import pandas as pd
import time
import multiprocessing
from timeit import default_timer as timer
from modulesT import getChildrenGoTerm

#%%

def power(x, n):

    # time.sleep(1)
    out = getChildrenGoTerm('test')
    return out

# values = ((2, 2), (4, 3), (5, 5))
query = ['a','b']
goTermsQuery = [['c','d','c','d','c','d','c','d'],['e','f','c','d','c','d','c','d','c','d']]
values = [[query[0],goTermsQuery],[query[1],goTermsQuery]]

start = timer()

if __name__ == '__main__':
    with multiprocessing.Pool(processes=2) as pool:
        results = pool.starmap(power, values)

end = timer()
print(end-start)
