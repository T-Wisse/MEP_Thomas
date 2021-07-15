---
title: "Thesis structure comments"
output: pdf_document
documentclass: article
classoption: onecolumn
pdf_document:
latex_engine: xelatex
toc: true
numberSections: true
highlight: tango
sectionsDepth: 3
chapters: True
figPrefix:
  - "Fig."
  - "Figs."
secPrefix:
  - "Section"
  - "Sections"
fontsize: 12pt
geometry: margin=0.5in
autoEqnLabels: true
cref: true
crossref: true
colorlinks: true
---

- Title: "*Exploring evolution using
Saturated Transposon
Analyis In Yeast*"

    - Comments:
        - The title does not convey the right picture of what you did in the thesis.
        - You should focus more on the GI part to understand a particular adaptive evolutionary pathway, like the one from $\Delta$ bem1 . 
        - An alternative could be: "*Genome-wise quantification of  genetic context influence on gene deletions in budding yeast* " **The case of the bem1-bem3 pair in the context of adaptive evolution**. 

1.  Introduction
    
    1. The $\Delta$bem1 adaptive evolutionary trajectory
    1. Genetic interactions in yeast
    1. State of the art (**of what??**) 
        - here you can also put , state of start approaches to measure GI in yeast. 
    1. SATAY as an alternative technique to compute GI
        - Description of the technique
        - why could it be a good candidate to quantify GI in yeast in a high-throughput manner. 
    1. Introduction to the research problem
        - Intro of what we still  dont know 
        - The relevance of tackling that question to gain understanding in how biological networks do change during evolution. 
    
    1. Research questions 
    1. Hypotheses
    1. Goals
        - Experimental goals
        - Modelling goals 



2. Materials & Methods
    1. Yeast strain construction
        - Construction of a bem1∆bem3∆ strain 
    1.  SATAY
        1. Library formation 
        2. DNA extraction 
        3. DNA sequencing
            - DNA Digestion
            - DNA ligation
            - PCR 
    1. Media
    1. OD measurements
        - nanodrop machine
    1. Population growth assays 
        - Biotek protocols

    1. Data analysis
        - Population growth assays
        - SATAY library complexity
        - Access to the software 
       

3. Model
    1. Current knowledge 
        - Interactors of bem1
        - Interactors of bem3
        - Interaction type between bem1 and bem3 , when bem3 is knocked down after a bem1 deletion. 
        
            - Change in fitness from population growth data from $\Delta$bem1 genotype to $\Delta$bem3$\Delta$bem1. 

        
        - Interaction type between bem1 and bem3 , when bem1 is knocked down after a bem3 deletion. 
            - Change in fitness from population growth data from $\Delta$bem3 genotype to $\Delta$bem3$\Delta$bem1. 

    1. What is predicted to change after bem1 deletion in WT and in $\Delta$bem3 background. 
        - focus on the genetic interactors from bem1 in both backgrounds
    1. What is predicted to change after bem3 deletion in WT and in $\Delta$bem1 background. 
        - focus on the genetic interactors from bem3 in both backgrounds

    1. Expanding on previous systems to know  which lethal deletions of $\Delta$bem3$\Delta$bem1 are not lethal in $\Delta$bem3
        - From Wessel results using ML on trying to predict gene essentiality using a different set of genes
        - Look at common observables that can be extrapolated to our system , like common go terms, common interactors , or functional enrichment of predicted genes in relation with the gene of interest. 

4. Results

    1. Experimental results
        - Complexity of the satay library 
        - Sequencing results 
        
    1. Modelling results 
        1. Prediction on what is the effect of bem1 deletion on a network that do not have bem3 , based on what is known in WT. 
        1.  Prediction of which type of genes are more prone to be SL in dbem1dbem3 background  that are not in the single knockouts backgrounds. 


5. Conclusion & Discussion

    - What do we learn after doing this research? 
    - Do we accomplish what we proposed to do at the beginning?
    - How our findings align with our initial hypotheses? 
    - Implications of the results for our new understanding on how biological networks do change during evolution. 

6. Future research 
    - What is still left and yet relevant to do? 
   
7. References 
8. Appendix
    - A Protocols 
    - B Primer sequences 
    - C Strains 
    - D Python code