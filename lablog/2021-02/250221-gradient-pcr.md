---
title: "Temperature gradient PCR"
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

# Title: Temperature gradient PCR
## Date
25-02-2021

## Method
### PCR
- primers OLIC22 and OLIC26: 100 uM 
- need 10uM: dilute 10x
- 60 ng DNA per reaction (1uL)
- vortex everything before pipetting. (carefully)
- buffer: 5x concentrated (phusion reaction buffer)
- dNTPs: 1uL

- adding a copy of yLL3a as control.
- dilute primer: 2uL in 18 uL milliQ

- prepare 50ml master mix:

| Amount         | Type          |
|----------------|---------------|
| 35.5 ul        | milliQ        |
| 10 ul          | buffer        |
| 1 ul           | dNTPs mix     |
| 1 ul           | primer1(10mM) |
| 1 ul           | primer2 (10mM)|
| 1 ul           | template DNA  |
| 0.5 ul         | Polymerase    |

- We will do 2 reactions, 1 for each template DNA -> make 100 ml
put 30uL (?) in the PCR tubes -> in the PCR machine

- Gel of PCR product shows no bands -> repeat with different annealing temperatures in the PCR

- Done a-h for 8 annealing temperatures:
65C, 64,3C, 63,1C, 59,0C, 57,3C, 56,0C, 55,0C
## Results
- (next day): DNA gel for the 8 samples: loading 10uL PCR product and 2uL loading dye. Ran at 120V for 20 min:
All lanes show bands
- Put an image of your gel
  
  `![](relative-path-to-the-image)`
  
## Conclusions