# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:45:10 2021

@author: Thomas
"""


#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# Type constraints should come early - before all mentions of the paths they constrain
query.add_constraint("goAnnotation.ontologyTerm", "GOTerm")
query.add_constraint("goAnnotation.ontologyTerm.parents", "GOTerm")

# The view specifies the output columns
query.add_view(
    "symbol", "featureType", "goAnnotation.ontologyTerm.name",
    "goAnnotation.ontologyTerm.namespace",
    "goAnnotation.ontologyTerm.parents.name"
)

# This query's custom sort order is specified below:
query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("goAnnotation.evidence.code.annotType", "=", "high-throughput", code="C")
query.add_constraint("goAnnotation.evidence.code.annotType", "=", "manually curated", code="D")
query.add_constraint("goAnnotation.qualifier", "!=", "NOT", code="I")
query.add_constraint("goAnnotation.qualifier", "IS NULL", code="J")
query.add_constraint("status", "=", "Active", code="E")
query.add_constraint("status", "IS NULL", code="F")
query.add_constraint("organism.name", "=", "Saccharomyces cerevisiae", code="B")
query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", "protein folding", code="A")

# Your custom constraint logic is specified with the code below:
query.set_logic("A and B and (C or D) and (E or F) and (I or J)")

for row in query.rows():
    print(row["symbol"], row["featureType"], row["goAnnotation.ontologyTerm.name"], \
        row["goAnnotation.ontologyTerm.namespace"], row["goAnnotation.ontologyTerm.parents.name"])