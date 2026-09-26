 # Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Bayesian estimation of deep-time diversification with PyRate](#1-Bayesian-estimation-of-deep-time-diversification-with-PyRate)
	- [1.1 P](#11-Preservation-model)
	- [1.2 Multivariate analyses](#12-Multivariateanalyses)
- [2 Phylogenetic comparative analyses](#2-Phylogenetic-comparative-analyses)
    - [2.1 Phylogenetic generalized linear regression](#21-Phylogeneticgeneralizedlinearregression)
	- [2.2 Analyses of discrete trait evolution with corHMM](#232-Analyses-of-discrete-trait-evolution-with-corHMM)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)" but can be generalized to other similar data. All of the presented scripts are written in R language (R Core Team, 2022).
	If you plan to use any of these scripts, please cite "XXX". </p>

## Overview

<p align="justify"> This repository contains html files for performing the following analyses:

**1**: 

**2**: Phylogenetic comparative analyses

<p align="justify"> All data used to perform each analysis are deposited in this repository </p>

## 1 Bayesian estimation of deep-time diversification with PyRate

`used directory (PyRate_scripts)`

<p align="justify">  In this first session, we will be using PyRate (Silvestro et al, 2014). PyRate is a program implemented in Python whose aim is to jointly estimate the preservation process, the tempo of origination and extinction of lineages based on their occurrences in the fossil record. Here, we will assume that the PyRate repository with its functions is at the root of the current working directory.</p>

### 1.1 Preservation model

`used directory (Preservation_Test)`

`used script (model_preservation_test.sh, model_drafting.r; run_preservation.sh)`

<p align="justify"> One of the main strengths of PyRate is its ability to account for the bias of the fossil record by estimating a preservation process and correcting the estimated age derived from raw occurrence data. Thus, choosing the best-fit preservation model for any PyRate analysis is critical. Fortunately, Silvestro et al. (2019) implemented a likelihood-based approach for preservation model selection. Yet, while this procedure is certainly useful, it is incomplete. Indeed, the first implementation allowed for model selection across HPP, NHPP, TPP and alternative versions of the TPP, with missing bins. However, bin removal occurred only once and was not recursive. Consequently, model selection is incomplete. Here, we corrected and enhanced this procedure by performing model selection on all PyRate replicates (here, 100). Furthermore, we allowed for recursive bin removal, meaning that the best fit TPP model could be a two-bin model, whereas the generating TPP model could be a five-bin model. Model selection is performed with pairwise comparisons of the AICc metrics across all replicates. </p>

### 1.2 Running PyRate

`used directory (PyRate_runs)`

`used script (PyRate_run.sh)`

<p align="justify"> The script provided in this section is rather simple, and runs a BDCS (birth-death with constrained shifts) analysis on 20,000,000 generations on the genus dataset, including singletons, with diversification shifts every 5 Myrs and integrating preservation shifts from the 1.1 section. Here, to be computationally efficient, we choose to parallelise our run on 20 CPUs </p> 

## 2 Phylogenetic comparative analyses

`used directory (Phylogenetic_comparative_analysis)`

<p align="justify"> In this section, we will perform several phylogenetic comparative analyses. </p>

### 2.1 Analyses of continuous trait evolution with OUwie and phylogenetic ANOVA

`used directories (Panova, OUwie)`

`used script (Phylogenetic_analysis_of_variance_(PANOVA).r; Phylogenetic_analysis_of_variance_(PANOVA-Replicated).r; PANOVA.sh; OUwie_consensus.r; OUwie_replicated.r; run_OUwie.sh)`

<p align="justify"> The second step in this section is to perform analyses of continuous trait evolution using both PANOVA and OUwie. The first step consists of running Phylogenetic analyses of variance to compare whether the mean difference between any number of groups significantly differs, even when considering phylogenetic relatedness. Here, to assess wether our results were robust regarding the analytical method employed, we implemented two phylogenetic ANOVA, the first with a null model process generated through simulation (sim-PANOVA; Revell, 2024), and a second with a null model process based on randomizing residuals in a permutation procedure (RRPP; Collyer & Adams, 2018). Both versions of these scripts are managed by the script "PANOVA.sh" which will perform PANOVA on all datasets and all trees (extant and fossil+extant). For OUwie, similarly to corHMM, both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Both of these scripts are managed by the "run_OUwie.sh" which essentially runs all these analyses on all trees (extant and fossil+extant) and traits (bioluminescence and habitat).  </p>


### 2.2 Analyses of discrete trait evolution with corHMM

`used directory (corHMM)`

`used script (corHMM_ASE_consensus.r; corHMM_ASE_replicated.r; run_corHMM.sh)`

<p align="justify"> The first step in this section is to perform analyses of discrete trait evolution using corHMM. Both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Both of these scripts are managed by the "run_corHMM.sh" which essentially runs all these analyses on all trees (extant and fossil+extant) and traits (bioluminescence and habitat).  </p>

### Reference

Brée, B., Condamine, F. L., & Guinot, G. (2022). Combining palaeontological and neontological data shows a delayed diversification burst of carcharhiniform sharks likely mediated by environmental change. Scientific Reports, 12(1), 21906.

Boyko, J. D., O’Meara, B. C., & Beaulieu, J. M. (2023). A novel method for jointly modeling the evolution of discrete and continuous traits. Evolution, 77(3), 836-851.

Collyer, M.L. & Adams, D.C. (2018) RRPP: an R package for fitting linear models to high-dimensional data using residual randomization. Methods in Ecology and Evolution, 9, 1772–1779.

Marion, A. F., Condamine, F. L., & Guinot, G. (2024). Sequential trait evolution did not drive deep-time diversification in sharks. Evolution, 78(8), 1405-1425.

Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ, 12, e16505.

R Core Team (2022). R: A language and environment for statistical computing. R Foundation for statistical computing, Vienna, Austria. URL https://www.R-project.org/.

Silvestro, D., Salamin, N., & Schnitzler, J. (2014). PyRate: a new program to estimate speciation and extinction rates from incomplete fossil data. Methods in Ecology and Evolution, 5(10), 1126-1131.

Silvestro, D., Salamin, N., Antonelli, A., & Meyer, X. (2019). Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Paleobiology, 45(4), 546-570.
