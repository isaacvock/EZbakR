# Welcome to EZbakR!

EZbakR is a highly flexible tool for analyses of nucleotide recoding RNA-seq datasets (NR-seq; e.g., [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582), [SLAM-seq](https://www.nature.com/articles/nmeth.4435), [TUC-seq](https://pubmed.ncbi.nlm.nih.gov/31768978/), etc.). See [our paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013179) for a discussion of the motivation behind EZbakR and its companion pipeline [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR), as well as validation of all of its novel functionality.

To install or update, run:

```
if (!require("roxygen2", quietly = TRUE))
    install.packages("roxygen2")
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("isaacvock/EZbakR")
```

At this point, changes will be made weekly, so updating frequently is highly recommended.

Documentation is here: [https://isaacvock.github.io/EZbakR/](https://isaacvock.github.io/EZbakR/)

### Want to learn more about how to interpret and analyze NR-seq data?

Check out some of the blogs I have written on the topic:

1. [An introduction on how to interpret NR-seq data](https://isaacvock.github.io/posts/post-with-code/). In this post, I use the task of building up an NR-seq simulator in R to help you develop an intuition for the important features of NR-seq data. Code is interactive and can be edited and re-executed in your browser without any additional setup (thanks to Quarto-live and webr).
2. [An introduction on how to analyze NR-seq data](https://isaacvock.github.io/posts/nrseq_analysis/). In this post, I motivate the concept of mixture modeling, the gold standard analysis strategy implemented in tools like EZbakR and GRAND-SLAM. I also show you how to implement the strategy yourself in R (code is interactive once again).
3. [An in-depth NR-seq review](https://isaacvock.github.io/posts/nrseq_deepdive/): Covers the history, data interpretation, data processing, analyzing processed data, designing experiments, and extensions of NR-seq that have been developed. This is taken from Chapter 1 of my [thesis](https://github.com/isaacvock/Thesis).

## Vignettes

Currently, the following functionalities have dedicated vignettes:

1. [Quickstart](https://isaacvock.github.io/EZbakR/articles/Quickstart.html): Takes you through the standard workflow, similar to bakR's one and only workflow.
2. [Estimating fractions](https://isaacvock.github.io/EZbakR/articles/EstimateFractions.html): Estimating the fraction of reads from each mutational population in your data. This is the nearly universal first step in all NR-seq analyses. This is done with EZbakR's [`EstimateFractions()`](https://isaacvock.github.io/EZbakR/reference/EstimateFractions.html) function.
3. [Estimating kinetics](https://isaacvock.github.io/EZbakR/articles/EstimateKinetics.html): Estimating kinetic parameters of synthesis and degradation in a standard NR-seq experiment. For standard, single label, NR-seq analyses, this is the next step in your analysis workflow after estimating the fraction of reads that are from labeled RNA. This is done with EZbakR's [`EstimateKinetics()`](https://isaacvock.github.io/EZbakR/reference/EstimateKinetics.html) function.
4. [Quality Control](https://isaacvock.github.io/EZbakR/articles/EZQC.html): Assessing the quality of your NR-seq data. This is done with EZbakR's [`EZQC()`](https://isaacvock.github.io/EZbakR/reference/EZQC.html) function.
5. [Comparative analyses](https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html): Fitting a flexible generalized linear model to your NR-seq data so as to perform comparative analyses of estimated kinetic parameters that complements differential expression analyses. This is done with EZbakR's [`AverageAndRegularize()`](https://isaacvock.github.io/EZbakR/reference/AverageAndRegularize.html) and [`CompareParameters()`](https://isaacvock.github.io/EZbakR/reference/CompareParameters.html) functions.
6. [Dynamical systems modeling](https://isaacvock.github.io/EZbakR/articles/EZDynamics.html). For analyses of subcellular fractionation and/or pre-mRNA processing dynamics. This is done with EZbakR's [`EZDynamics()`](https://isaacvock.github.io/EZbakR/reference/EZDynamics.html) function.
7. [Navigating EZbakR output](https://isaacvock.github.io/EZbakR/articles/EZget.html). Conveniently fetching data from EZbakR analyses. This is done with EZbakR's [`EZget()`](https://isaacvock.github.io/EZbakR/reference/EZget.html) function.

Other implemented functionality that may be of interest includes:

1. Providing fractions or kinetic parameter estimates as input. The former works similarly to how it did in bakR, and is implemented via the [`EZbakRFractions()`](https://isaacvock.github.io/EZbakR/reference/EZbakRFractions.html) function. THe latter is unique to EZbakR and is implemented via the [`EZbakRKinetics()`](https://isaacvock.github.io/EZbakR/reference/EZbakRKinetics.html) function.
2. Simulating NR-seq data. There are a number of simulation functions implemented in EZbakR. [`EZSimulate()`](https://isaacvock.github.io/EZbakR/reference/EZSimulate.html) is a convenient wrapper to several of these.

## Update (3/15/2025): Analyses of transcript isoforms

We recently (3/14/2025) put out [a preprint](https://www.biorxiv.org/content/10.1101/2025.03.12.642874v1) describing a method by which to analyze the kinetics of individual transcript isoforms using short read NR-seq data from total RNA. While this strategy is touched on a little bit in one of the EZbakR vignettes ([this one](https://isaacvock.github.io/EZbakR/articles/EstimateFractions.html#isoform-deconvolution)), I have also developed a full fastq-to-volcano plot walkthrough using real downsampled fastq files from that preprint so you can see how every step of the fastq2EZbakR and EZbakR pipeline needs to be configured/run for these analyses. The tutorial is [here](https://isaacvock.github.io/Isoform_Tutorial_Docs/), and the data used in that tutorial is [here](https://github.com/isaacvock/Isoform_Analysis_Tutorial). Over the next couple weeks I will be adding some extra details/analyses to this tutorial, but in its current form (as of 3/15/2025), all of the basics of performing isoform-level analyses are covered there. It also acts as a hand-on tutorial for all of the EZbakR-suite and can thus useful to checkout and try out even if you aren't interested in this particular analysis strategy.

## What's new?

EZbakR represents a complete rewrite of [bakR](https://github.com/simonlabcode/bakR). Improvements implemented in EZbakR include:

1. Modular function design that facilitates using EZbakR with any kind of NR-seq data, regardless of the experimental design or data details.
2. Extended mixture modeling capabilities. Includes:
    * Support for multi-label analyses.
    * Hierarchical mutation rate estimation strategy to allow for feature-specific mutation rates.
    * More efficient and accurate uncertainty quantification.
3. Additional kinetic parameter estimation strategies:
    * Non-steady-state analyses as introduced in [Narain et al., 2021](https://www.sciencedirect.com/science/article/pii/S1097276521004962).
    * Short-feed analyses that assume negligible degradation of existing RNA.
    * Synthesis rate estimation is implemented as a part of all strategies.
4. Improved uncertainty propogation so as to achieve performance of bakR's slower implementations (Hybrid and MCMC) with a strategy as efficient as bakR's most efficent implementation (MLE).
5. Removal of Stan dependencies. I love Stan, but having it as an R package dependency makes installation and maintenace more difficult. 
6. Optional [Apache Arrow](https://arrow.apache.org/) backend to help with analyses of larger-than-RAM datasets
7. Linear model-based averaging of replicate data to support more complex experimental designs and maximally flexible comparative analyses. 
8. Greater flexibility in terms of the input data structure. Namely, multiple different features can be specified in your input cB table, and multiple different experimental details can be included in your input metadf table.
9. A novel transcript isoform deconvolution strategy that allows for isoform-specific kinetic parameter estimation.
10. Generalized linear dynamical systems modeling of NR-seq data. Supports analyses of subcellular fractionation NR-seq extensions, such as those described [here](https://www.cell.com/molecular-cell/fulltext/S1097-2765(24)00511-2#:~:text=Thus%2C%20RNA%20flow%20impacts%20cell,processing%2C%20including%20splicing%20and%20polyadenylation.), [here](https://www.biorxiv.org/content/10.1101/2024.03.11.584215v1), and [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012059). Also supports analyses of pre-mRNA processing dynamics.

In the near future, EZbakR will support anything bakR can do that isn't currently implemented (Namely `DissectMechanisms()` and various visualization functions). There are also a number of exciting developments on the horizon, so stay tuned!

## What is NR-seq?

NR-seq refers to a class of methods that combine RNA-seq, metabolic labeling, and unique metabolic label recoding chemistries. These methods were originally developed to dissect the kinetics of RNA synthesis and degradation. Excitingly though, a treasure trove of extensions of the original methods have been created over the years. To-date, nucleotide recoding has been combined with the likes of [TT-seq](https://www.nature.com/articles/nmeth.4582), [Start-seq](https://www.sciencedirect.com/science/article/pii/S1097276521006869?via%3Dihub), [Ribo-seq](https://www.nature.com/articles/s41592-021-01250-z), [scRNA-seq](https://www.nature.com/articles/s41586-019-1369-y) (other examples of this [here](https://www.nature.com/articles/s41592-020-0935-4), [here](https://www.nature.com/articles/s41587-020-0480-9), and [here](https://www.biorxiv.org/content/10.1101/2023.07.06.547989v1)), [Perturb-seq](https://www.nature.com/articles/s41587-023-01948-9), [long-read sequencing](https://www.biorxiv.org/content/10.1101/2020.05.01.073296v1), and [subcellular fractionation](https://www.biorxiv.org/content/10.1101/2022.08.21.504696v1.full). In addition, while the original methods used 4-thiouridine (s<sup>4</sup>U), the same chemistry has been found to work with [6-thioguanosine](https://pubs.acs.org/doi/full/10.1021/jacs.8b08554) (s<sup>6</sup>G), opening the door to dual-labeling experimental designs (e.g., [TILAC](https://academic.oup.com/nar/article/50/19/e110/6677324)). EZbakR and its companion pipeline [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR) aim to provide an integrated and flexible framework to support this exciting class of methods. 
