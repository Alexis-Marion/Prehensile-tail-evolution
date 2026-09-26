# Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Preliminary analyses](#1-Preliminaryanalyses)
	- [1.1 Raw data analyses](#11-Rawdataanalyses)
	- [1.2 Multivariate analyses](#12-Multivariateanalyses)
- [2 Phylogenetic comparative analyses](#2-Phylogenetic-comparative-analyses)
    - [2.1 Phylogenetic generalized linear regression](#21-Phylogeneticgeneralizedlinearregression)
	- [2.2 Analyses of discrete trait evolution with corHMM](#22-Analyses-of-discrete-trait-evolution-with-corHMM)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)" but can be generalized to other similar data. All of the presented scripts are written in R language (R Core Team, 2022). Most of the functions employed here are featured in the geomorph package (Adams& Otárola‐Castillo, 2013), phytools 2.0 (Revell, 2024) and CorHMM (Boyko & Beaulieu, 2021). If you plan to use any of these scripts, please cite "XXX". </p>

## Overview

<p align="justify"> This repository contains html files for performing the following analyses:

**1**: Multivariate analyses

**2**: Phylogenetic comparative analyses

<p align="justify"> All data used to perform each analysis are deposited in this repository </p>

## 1 Preliminary analyses

`used directory (Preliminary analyses)`

<p align="justify">  In this first session, we will be mainly focusing on multivariate analyses of morphological data.</p>

### 1.1 Raw data analyses

`used script (Raw_data.ipynb)`

<p align="justify"> The Purpose of this first script is to load and clean raw data obtained through landmarking procedures. In this script, based on the tail length and landmarks, two main metrics will be computed : the Transverse ProcEss Index (TPEI) and the Robusticity Index (RI). Each of these metrics will be computed for all vertebrae, distal vertebrae, transitional vertebrae, proximal vertebrae and the last 25% vertebrae remaining (aka Last Quarter). Species (or lineages) not represented in the phylogenetic tree are removed from the dataset. All metrics are computed for each species and merged with ecological information data in a synthetic dataset.</p>

### 1.2 Multivariate analyses

`used script (Multivariate analyses & simple phylogenetic regression.ipynb)`

<p align="justify"> The Purpose of this first script is to provide a collection of multivariate analyses and some basic phylogenetic comparative analyses of the data computed in subsection 1.1. In this script, one can perform Principal Component Analysis (PCA), Linear Discriminent Analyses to examine whether certain tail morphologies clustered together relative to their species ecology. Phylogenetic comparative analyses, such as phylogenetic regression aim at examining wether certain morphologies associated with tail prehensility drift from isometry expectation Lastly, convergence indices are computed. </p>

<p align="justify"> </p> 

## 2 Phylogenetic comparative analyses

`used directory (Phylogenetic_comparative_analyses)`

<p align="justify"> In this section, we will perform several phylogenetic comparative analyses. </p>

### 2.1 Phylogenetic generalised linear regression and Phylogenetic ANCOVA

`used script (Pgls Pancova.ipynb)`

<p align="justify"> The purpose of this script is to provide a formal phylogenetic-informed assessment of caudal morphological differences between ecological categories. To do so, generalised regression per ecological categories are estimated, and phylogenetic analyses of covariances (between index value and size) are computed. </p>

### 2.2 Analyses of discrete trait evolution with corHMM

`used script (Ancestral_state_estimation_PT.ipynb, Ancestral_state_estimation_PT_replicated.ipynb)`

<p align="justify"> In this section analyses of discrete trait evolution are performed using corHMM. Both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Based on theses results, estimates of ancestral trait are computed, and displayed directly onto the tree of interest. </p>

### Reference

Adams, D. C., & Otárola‐Castillo, E. (2013). geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in ecology and evolution, 4(4), 393-399.

Boyko, J. D., & Beaulieu, J. M. (2021). Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution, 12(3), 468-478.

Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ, 12, e16505.

R Core Team (2022). R: A language and environment for statistical computing. R Foundation for statistical computing, Vienna, Austria. URL https://www.R-project.org/.
