 # Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Bayesian estimation of deep-time diversification with PyRate](#1-Bayesian-estimation-of-deep-time-diversification-with-PyRate)
	- [1.1 Preservation model](#11-Preservation-model)
	- [1.2 Running PyRate](#12-Running-Pyrate)
	- [1.3 Assess convergence](#13-Assess-convergence)
  	- [1.4 Plotting PyRate results](#14-Plotting-PyRate-results)
- [2 Analyses of the fossil record](#2-Analyses-of-the-fossil-record)
 	- [2.1 Selecting appropriate PyRate models](#21-Selecting-appropriate-PyRate-models)
	- [2.2 Extracting time for speciation and extinction](#22-Extracting-time-for-speciation-and-extinction)
	- [2.3 Estimating lineage through time per habitat](#23-Estimating-lineage-through-time-per-habitat)
  	- [2.4 Estimating tempo of origination](#24-Estimating-tempo-of-origination)
	- [2.5 Grafting fossils](#25-Grafting-fossils)
- [3 Phylogenetic comparative analyses](#3-Phylogenetic-comparative-analyses)
	- [3.1 Analyses of discrete trait evolution with corHMM](#31-Analyses-of-discrete-trait-evoplution-with-corHMM)
 	- [3.2 Analyses of continuous trait evolution with OUwie and phylogenetic ANOVA](#32-Analyses-of-continuous-trait-evolution-with-OUwie-and-phylogenetic-ANOVA)
  	- [3.3 Joint estimation of discrete and continuous trait evolution with hOUwie](#32-Joint-estimation-of-discrete-and-continuous-trait-evolution-with-hOUwie)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Caudal vertebral convergence in relation to tail prehensility in Murinae (Rodentia)" but can be generalized to other similar data. All of the presented scripts are written in R language (R Core Team, 2022).
	If you plan to use any of these scripts, please cite "XXX". </p>

## Overview

<p align="justify"> This repository contains html files for performing the following analyses:

**1**: 

**2**: Analyses of the fossil record

**3**: Phylogenetic comparative analyses

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

### 1.3 Assess convergence

`used directory (PyRate_convergence)`

`used script (assess_run_convergence.py; plot_ess.r; run_convergence.sh)`

<p align="justify"> PyRate is a Bayesian program, thus, we will consider that a PyRate run is finished when it achieves convergence. A popular metric to evaluate convergence is the ESS (effective sample size), and it is generally considered that when its number is above 200, convergence is achieved. Thus, we assessed convergence on all runs using the "assess_run_convergence.py" script. Furthermore, this script gives additional useful metrics on the run, such as the origination age, extinction age (if including solely extinct taxa), and the proportion of Ts and Te, with ESS above 200. We also provided a graphical output that can be executed using "plot_ess.r". These two scripts can be run sequentially using the "run_convergence.sh" script.

### 1.4 Plotting PyRate results

`used directory (PyRate_plotting)`

`used script (1-extract_param_from_PyRate_outputs.r; 2-plotting_facilities.r; Plot_rates.r; Q_rate.sh; ltt_creator.sh; master_script_plotting.sh; parse_Q_rates.py; run_plotting.sh)`

<p align="justify"> In this last section, we provided plotting scripts to display graphically each PyRate output. These scripts will take as input the output directory of a regular (BDCS or RJMCMC) PyRate run. They will represent the RTT (origination and extinction), the diversification RTT, the LTT and the QTT. All the aforementioned scripts are managed using the "master_script_plotting.sh" script, which is the only one that is meant to be run by the user. </p>

## 2 Analyses of the fossil record

`used directory (Analyses_of_the_fossil_record)`

<p align="justify">  In this section session, we will be using the output from PyRate and molecular phylogenies to estimate the tempo of origination of several clades, ecologies and traits.</p>

### 2.1 Selecting appropriate PyRate models

`used directory (Model_selection)`

`used script (Model_selection.r; run_Model_selection.sh)`

<p align="justify"> In this section, we propose a quantitative comparison between each PyRate model couple (e.g. a model with full preservation time bins, and another with fewer preservation shift frames). Specifically, since some alternative models do not possess the same structure, statistical comparison between each model could be inappropriate. Instead, we could still compare the estimated values (in our case Ts, Te, origination and extinction rate) and evaluate whether they significantly differ when comparing two alternative models. Here we choose to compare each PyRate replicate separately (thus each comparison is performed 20 times), and consider performing pairwise comparisons of Ts, Te, origination and extinction rate values between each two alternative models. </p>

### 2.2 Extracting time for speciation and extinction

`used directory (Ts_Te_maker)`

`used script (Ts_Te_combined.r)`

<p align="justify"> The script used in this section permits to create additional Ts_Te files that can be later used in PyRate. The combination method used here is that of Brée et al. (2022) and consists of extracting the estimated Time for Speciation (Ts) and Time for Extinction (Te) estimated with PyRate, and merging them with the estimated extant taxa age from a dated phylogeny. Here, we did so by extracting the branch length  from 100 posterior trees (Marion et al., 2024) of all extant species that are not sampled in the fossil record and merging them with the Ts Te estimated with PyRate.  </p>

### 2.3 Estimating lineage through time per habitat

`used directory (Lineage_through_time)`

`used script (Lineage through time by habitat (Squaliformes).r; Lineage through time by habitat (Squaliformes)-Randomised.r; utility_function_ltt_plot.r)`

<p align="justify"> Using the Ts/Te output from section 2.2 and by using palaeoenvironmental data, we estimated the lineage through time for each habitat. The first script "Lineage through time by habitat (Squaliformes).r" plot the lineage through time per habitat for each dataset, whereas the second script "(Lineage through time by habitat (Squaliformes)-Randomised.r" randomly shuffle the habitat identity of each lineage while keeping the lineage through time pattern unaltered. Both these scripts use the "utility_function_ltt_plot.r" script.  </p>

### 2.4 Estimating tempo of origination

`used directory (Origination_tempo)`

`used script (Plot_estimated_ages.r)`

<p align="justify"> Using the Ts/Te output from section 2.2 and the posterior tree distribution, this script estimates the tempo of origination for each major squaliform clade, using both palaeontological and neontological evidence. Likewise, this script also estimates the tempo of origination of each habitat preference across all Squaliformes families, using the Ts/Te values as input. Lastly, by combining palaeontological and neontological evidence, we estimated the tempo of origination of each putative bioluminescent clade. </p>

### 2.5 Grafting phylogenies

`used directory (Grafting_Trees)`

`used script (Grafting_trees.r; function_grafting.r)`

<p align="justify"> In this last section, we provide tools to graft fossil taxa onto extant phylogenies. Placing a fossil into an already resolved extant time-calibrated tree represents a major challenge. To account for the potential bias possibly affecting the placement of fossils, we relied on previously published phylogenetic hypotheses. Furthermore, to account for temporal uncertainty, we randomly sampled the origination age for each extinct lineage from a uniform distribution bounded with the estimated PyRate origination age as the upper limit, and the stem age of the sister clade of the fossil as the lower limit (i.e. the node preceding the sister-clade mentioned above). These scripts allow us to perform this procedure on the consensus and posterior distribution, grafting up to twelve fossils. </p>

## 3 Phylogenetic comparative analyses

`used directory (Phylogenetic_comparative_analysis)`

<p align="justify"> In this section, we will perform several phylogenetic comparative analyses. </p>

### 3.1 Analyses of discrete trait evolution with corHMM

`used directory (corHMM)`

`used script (corHMM_ASE_consensus.r; corHMM_ASE_replicated.r; run_corHMM.sh)`

<p align="justify"> The first step in this section is to perform analyses of discrete trait evolution using corHMM. Both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Both of these scripts are managed by the "run_corHMM.sh" which essentially runs all these analyses on all trees (extant and fossil+extant) and traits (bioluminescence and habitat).  </p>

### 3.2 Analyses of continuous trait evolution with OUwie and phylogenetic ANOVA

`used directories (Panova, OUwie)`

`used script (Phylogenetic_analysis_of_variance_(PANOVA).r; Phylogenetic_analysis_of_variance_(PANOVA-Replicated).r; PANOVA.sh; OUwie_consensus.r; OUwie_replicated.r; run_OUwie.sh)`

<p align="justify"> The second step in this section is to perform analyses of continuous trait evolution using both PANOVA and OUwie. The first step consists of running Phylogenetic analyses of variance to compare whether the mean difference between any number of groups significantly differs, even when considering phylogenetic relatedness. Here, to assess wether our results were robust regarding the analytical method employed, we implemented two phylogenetic ANOVA, the first with a null model process generated through simulation (sim-PANOVA; Revell, 2024), and a second with a null model process based on randomizing residuals in a permutation procedure (RRPP; Collyer & Adams, 2018). Both versions of these scripts are managed by the script "PANOVA.sh" which will perform PANOVA on all datasets and all trees (extant and fossil+extant). For OUwie, similarly to corHMM, both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Both of these scripts are managed by the "run_OUwie.sh" which essentially runs all these analyses on all trees (extant and fossil+extant) and traits (bioluminescence and habitat).  </p>


### 3.3 Joint estimations of discrete and continuous trait evolution with hOUwie

`used directory (hOUwie)`

`used script (hOUwie_consensus.r; hOUwie_consensus_bioluminescence_Replicated.r; run_hOUwie.sh)`

<p align="justify"> For this last section, we will perform joint estimation of discrete and continuous trait evolution with hOUwie (Boyko & Beaulieu, 2023). Like corHMM and OUwie, both versions of this script (consensus vs replicated) designate the consensus tree and the posterior distribution, respectively. Both of these scripts are managed by the "run_hOUwie.sh", which runs all these analyses on traits (bioluminescence and habitat). </p>


### Reference

Brée, B., Condamine, F. L., & Guinot, G. (2022). Combining palaeontological and neontological data shows a delayed diversification burst of carcharhiniform sharks likely mediated by environmental change. Scientific Reports, 12(1), 21906.

Boyko, J. D., O’Meara, B. C., & Beaulieu, J. M. (2023). A novel method for jointly modeling the evolution of discrete and continuous traits. Evolution, 77(3), 836-851.

Collyer, M.L. & Adams, D.C. (2018) RRPP: an R package for fitting linear models to high-dimensional data using residual randomization. Methods in Ecology and Evolution, 9, 1772–1779.

Marion, A. F., Condamine, F. L., & Guinot, G. (2024). Sequential trait evolution did not drive deep-time diversification in sharks. Evolution, 78(8), 1405-1425.

Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ, 12, e16505.

R Core Team (2022). R: A language and environment for statistical computing. R Foundation for statistical computing, Vienna, Austria. URL https://www.R-project.org/.

Silvestro, D., Salamin, N., & Schnitzler, J. (2014). PyRate: a new program to estimate speciation and extinction rates from incomplete fossil data. Methods in Ecology and Evolution, 5(10), 1126-1131.

Silvestro, D., Salamin, N., Antonelli, A., & Meyer, X. (2019). Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Paleobiology, 45(4), 546-570.
