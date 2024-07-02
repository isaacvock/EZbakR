# Welcome to the alpha release of EZbakR!

EZbakR is a one-stop-shop for analyses of nucleotide recoding RNA-seq datasets(NR-seq). NR-seq refers to a class of methods (e.g., [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582), [SLAM-seq](https://www.nature.com/articles/nmeth.4435), [TUC-seq](https://pubmed.ncbi.nlm.nih.gov/31768978/), etc.) that combine RNA-seq, metabolic labeling, and unique metabolic label recoding chemistries. These methods were originally developed to dissect the kinetics of RNA synthesis and degradation. Excitingly though, a treasure trove of extensions of the original methods have been created over the years. To-date, nucleotide recoding has been combined with the likes of [TT-seq](https://www.nature.com/articles/nmeth.4582), [Start-seq](https://www.sciencedirect.com/science/article/pii/S1097276521006869?via%3Dihub), [Ribo-seq](https://www.nature.com/articles/s41592-021-01250-z), [scRNA-seq](https://www.nature.com/articles/s41586-019-1369-y) (other examples of this [here](https://www.nature.com/articles/s41592-020-0935-4), [here](https://www.nature.com/articles/s41587-020-0480-9), and [here](https://www.biorxiv.org/content/10.1101/2023.07.06.547989v1)), [Perturb-seq](https://www.nature.com/articles/s41587-023-01948-9), [long-read sequencing](https://www.biorxiv.org/content/10.1101/2020.05.01.073296v1), and [subcellular fractionation](https://www.biorxiv.org/content/10.1101/2022.08.21.504696v1.full). In addition, while the original methods used 4-thiouridine (s<sup>4</sup>U), the same chemistry has been found to work with [6-thioguanosine](https://pubs.acs.org/doi/full/10.1021/jacs.8b08554) (s<sup>6</sup>G), opening the door to dual-labeling experimental designs (e.g., [TILAC](https://academic.oup.com/nar/article/50/19/e110/6677324)).

To install or update, run:

```
library(devtools)
devtools::install_github("isaacvock/EZbakR")
```

At this point, changes will be made weekly, so updating frequently is highly suggested.

Documentation is here: https://isaacvock.github.io/EZbakR/

## What's new?

EZbakR represents a complete rewrite of bakR. Improvements implemented in EZbakR include:

1. Modular function design that facilitates using EZbakR with any kind of NR-seq data, regardless of the experimental design or data details.
2. Extended mixture modeling capabilities. Includes:
    * Support for multi-label analyses.
    * Hierarchical mutation rate estimation strategy to allow for feature-specific mutation rates.
    * More efficient and accurate uncertainty quantification
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

In the near future, EZbakR will support:
1. Alternative overdisperse mixture modeling strategies.
2. Generalized dynamical systems modeling of NR-seq data.
3. Anything bakR can do that isn't currently implemented (Namely `DissectMechanisms()`, `QC_checks()`, and various visualization functions).
