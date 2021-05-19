---
title: "random"
output: pdf_document
documentclass: article
classoption: twocolumn
pdf_document:
latex_engine: pdflatex
toc: true
lof: true
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

# PCR of genomic DNA from yLL117a  

## Date 
05/03/21 - 09/03/21

## Objective
Obtain a sufficient concentration of DNA containing BEM1d

## Method
Incubation of yeast cells (yLL117a) in liquid YPD for several days until a dense culture has formed. Incubating 2*10ml and negative control.
Genomic DNA extraction of yLL17a popultion using kit (GeneMATRIX gram plus & yeast genomic DNA purification kit). Using ~18 ml worth of pellet over 10 epis. Using 30uL Elution buffer per epi (80C). Combinined all DNA into 1 epi.
Checked extracted DNA concentration with nanodrop.

Performed PCR on extracted yLL117a DNA and previously extracted yLL3a DNA using primers OLIC54 and OLIC55 on both. These primers flank the gene BEM1.
Used volume: 50 uL per PCR tube. 8 tubes per DNA using a temperature gradient from 55C to 65C.
Large gel (110V for 30min) to check presence, length and concentration of  in PCR product. Loaded 10uL per lane, except gDNA (5uL)

## Results
Nanodrop measurement of gDNA yLL117a:
60.4 ng/uL. 260/280: 1.88, 260/230: 1.16

![a longer caption for testing](../Images/210309-PCR_yLL117_olic54_oli55_annotated.png)
Caption:
Figure 1: DNA gel with PCR products using primers OLIC54 and OLIC55 (flanking BEM1). DNA is present in all lanes (Measure concentration with biorad...). yLL117a PCR bands show a length of ~2.5kb. yLL3a PCR bands (control) show a length between 2.5kb and 3kb. 

## Conclusion
Based on the nanodrop measurement and gel, the DNA extraction of yLL117a was succesful. A 260/280 of 1.88 suggests the DNA is not pure (~1.8 is pure DNA), but also contains RNAs (2.0 is pure RNA). 260/230 of 1.16 (OR 2.16??) would also indicate contamination (2.0-2.2 expected for pure Nucleic acids)
The gel shows PCR product for both yLL3a and yLL117a as expected and for all annealing temperatures. yLL3a contains BEM1 and with primers OLIC54 and 55 should have a length of 2.7kb. yLL117a contains a construct in place of BEM1 and and with primers OLIC54 and 55 should have a length of 2.5 kb. These expected lengths correspond with observed lengths. This suggests we have sucessfully obtained PCR product of a BEM1delete construct from yLL117a.

