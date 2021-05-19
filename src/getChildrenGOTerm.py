# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 15:33:29 2021

@author: Thomas
"""

#!/usr/bin/env python
from __future__ import print_function
from intermine.webservice import Service
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
query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", "establishment or maintenance of cell polarity", code="A")
query.set_logic("A and B and (C or D) and (E or F) and (I or J)")

for row in query.rows():
    print(row["goAnnotation.ontologyTerm.name"], row["goAnnotation.ontologyTerm.identifier"])