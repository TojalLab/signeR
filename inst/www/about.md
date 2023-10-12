---
title: "About"
output: html_document
---

The signeRFlow is an interactive tool for exploring and analyzing signeR methods and TCGA signature data.

# About signeRFlow

signeRFlow is a suite of algorithms and datasets organized to allow the exploration of mutational signatures and exposures to related mutational processes. 

Available tools permit the analysis of user data or the exploration of public (TCGA) data.

Exposure data may be used to cluster samples,  and their relation to other sample data, such as clinical or survival data, can also be explored.

## How to cite us

Please cite the [signeR](https://academic.oup.com/bioinformatics/article/33/1/8/2525683?login=false) paper (Pubmed: [27591080](https://pubmed.ncbi.nlm.nih.gov/27591080/)) when you use the signeRFlow app.
```md
Rafael A Rosales, Rodrigo D Drummond, Renan Valieris, Emmanuel Dias-Neto, Israel T da Silva, signeR: an empirical Bayesian approach to mutational signature discovery, Bioinformatics, Volume 33, Issue 1, 1 January 2017, Pages 8–16, https://doi.org/10.1093/bioinformatics/btw572
```

## signeRFlow 

The main feature of the signeRFlow shiny app is the signeR modules, which provide access to signeR analysis to explore and visualize results interactively. Each module presents parameters and information, with views to all the available plots in the signeR package.

 - **signeR *de novo*:** This module provides access to signeR de novo
                analysis to find signatures in your data,
                estimating both signatures and related exposures.
 - **signer fitting**: This module provides access to signeR fitting
                analysis to find exposures to known signatures in your data, which can be uploaded or chosen from thethe  Cosmic database.

## TCGA Explorer

Another feature of the signeRFlow shiny app is the TCGA Explorer, which provides access to the results of signeR applications to TCGA datasets (33 cancer types). Here, we used the MC3 available mutations to run signeR across 33 different cancer types.

The [MC3](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub) (*Multi-Center Mutation Calling in Multiple Cancers*) project is an effort to generate a comprehensive encyclopedia of somatic mutation calls for the TCGA data to enable robust cross-tumor-type analyses. 

In this module, you can explore the results of signeR analysis using *de novo* and fitting modules.
