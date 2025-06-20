% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulations.R
\name{EZSimulate}
\alias{EZSimulate}
\title{Simulate NR-seq data for multiple replicates of multiple biological conditions}
\usage{
EZSimulate(
  nfeatures,
  mode = c("standard", "dynamics"),
  ntreatments = ifelse(mode == "standard", 2, 1),
  nreps = 3,
  nctlreps = 1,
  metadf = NULL,
  mean_formula = NULL,
  param_details = NULL,
  seqdepth = nfeatures * 2500,
  label_time = 2,
  pnew = 0.05,
  pold = 0.001,
  readlength = 200,
  Ucont = 0.25,
  feature_prefix = "Gene",
  dispslope = 5,
  dispint = 0.01,
  logkdegsdtrend_slope = -0.3,
  logkdegsdtrend_intercept = -2.25,
  logksynsdtrend_slope = -0.3,
  logksynsdtrend_intercept = -2.25,
  logkdeg_mean = -1.9,
  logkdeg_sd = 0.7,
  logksyn_mean = 2.3,
  logksyn_sd = 0.7,
  logkdeg_diff_avg = 0,
  logksyn_diff_avg = 0,
  logkdeg_diff_sd = 0.5,
  logksyn_diff_sd = 0.5,
  pdiff_kd = 0.1,
  pdiff_ks = 0,
  pdiff_both = 0,
  pdo = 0,
  dynamics_preset = c("preRNA", "nuc2cyto", "preRNAwithPdeg", "nuc2cytowithNdeg",
    "subtlseq", "nuc2cytowithpreRNA"),
  unassigned_name = "__no_feature",
  dispersion = 1000,
  lfn_sd = 0.2,
  treatment_effects = NULL,
  effect_avg_default = 0,
  effect_sd_default = 0.5,
  fraction_affected_default = 0.5,
  log_means = NULL,
  log_sds = NULL
)
}
\arguments{
\item{nfeatures}{Number of "features" (e.g., genes) for which to simulated data.}

\item{mode}{Currently, EZSimulate can simulate in two modes: "standard" and "dynamics".
The former is the default and involves simulating multiple conditions of standard NR-seq data.
"dynamics" calls \code{SimulateDynamics()} under the hood to simulate a dynamical systems
model of your choice. Most of the additional parameters do not apply if mode == "dynamics",
except for those from dynamics_preset and on.}

\item{ntreatments}{Number of distinct treatments to simulate. This parameter is
only relevant if \code{metadf} is not provided.}

\item{nreps}{Number of replicates of each treatment to simulate. This parameter is
only relevant if \code{metadf} is not provided}

\item{nctlreps}{Number of -s4U replicates of each treatment to simulate. This parameter
is only relevant if \code{metadf} is not provided.}

\item{metadf}{A data frame with the following columns:
\itemize{
\item sample: Names given to samples to simulate.
\item \if{html}{\out{<details>}}: Any number of columns with any names (not taken by other metadf columns)
storing factors by which the samples can be stratified. These can be referenced
in \code{mean_formula}, described below.
}
These parameters (described more below) can also be included in metadf to specify sample-specific simulation
parameter:
\itemize{
\item seqdepth
\item label_time
\item pnew
\item pold
\item readlength
\item Ucont
}}

\item{mean_formula}{A formula object that specifies the linear model used to
relate the factors in the \if{html}{\out{<details>}} columns of \code{metadf} to average log(kdegs) and
log(ksyns) in each sample.}

\item{param_details}{A data frame with one row for each column of the design matrix
obtained from \code{model.matrix(mean_formula, metadf)} that describes how to simulate
the linear model parameters. The columns of this data frame are:
\itemize{
\item param: Name of linear model parameter as it appears in the column names of the
design matrix from \code{model.matrix(mean_formula, metadf)}.
\item reference: Boolean; TRUE if you want to treat that parameter as a "reference". This
means that all other parameter values that aren't global parameters are set equal to this
unless otherwise determined (see \verb{pdiff_*} parameters for how it is determined if a parameter
will differ from the reference).
\item global: Boolean; TRUE if you want to treat that parameter as a global parameter. This means
that a single value is used for all features.
\item logkdeg_mean: If parameter is the reference, then its value for the log(kdeg) linear model
will be drawn from a normal distribution with this mean. If it is a global parameter, then this
value will be used. If it is neither of these, then its value in the log(kdeg) linear model will
either be the reference (if there is no difference between this condition's value and the reference)
or the reference's value + a normally distributed random variable centered on this value.
\item logkdeg_sd: sd used for draws from normal distribution as described for \code{logkdeg_mean}.
\item logksyn_mean: Same as \code{logkdeg_mean} but for log(ksyn) linear model.
\item logksyn_sd: Same as \code{logkdeg_sd} but for log(kdeg) linear model.
\item pdiff_ks: Proportion of features whose value of this parameter in the log(ksyn) linear model
will differ from the reference's. Should be a number between 0 and 1, inclusive. For example, if
\code{pdiff_ks} is 0.1, then for 10\% of features, this parameter will equal the reference parameter +
a normally distributed random variable with mean \code{logksyn_mean} and sd \code{logksyn_sd}. For the other
90\% of features, this parameter will equal the reference.
\item pdiff_kd: Same as \code{pdiff_ks} but for log(kdeg) linear model.
\item pdiff_both: Proportion of features whose value for this parameter in BOTH the
log(kdeg) and log(ksyn) linear models will differ from the reference. Value must be
between 0 and min(c(pdiff_kd, pdiff_ks)) in that row.
}
If param_details is not specified by the user, the first column of the design matrix
is assumed to represent the reference parameter, all parameters are assumed to be
non-global, logkdeg_mean and logksyn_mean are set to the equivalently named parameter values
described below for the reference and \code{logkdeg_diff_avg} and \code{logksyn_diff_avg} for all other parameters,
logkdeg_sd and logksyn_sd are set to the equivalently named parameter values
described below for the reference and \code{logkdeg_diff_sd} and \code{logksyn_diff_sd} for all other parameters,
and pdiff_kd, pdiff_ks, and pdiff_both are all set to the equivalently named parameter values.}

\item{seqdepth}{Total number of reads in each sample.}

\item{label_time}{Length of s^4^U feed to simulate.}

\item{pnew}{Probability that a T is mutated to a C if a read is new.}

\item{pold}{Probability that a T is mutated to a C if a read is old.}

\item{readlength}{Length of simulated reads. In this simple simulation, all reads
are simulated as being exactly this length.}

\item{Ucont}{Probability that a nucleotide in a simulated read is a U.}

\item{feature_prefix}{Name given to the i-th feature is \code{paste0(feature_prefix, i)}. Shows up in the
\code{feature} column of the output simulated data table.}

\item{dispslope}{Negative binomial dispersion parameter "slope" with respect to read counts. See
DESeq2 paper for dispersion model used.}

\item{dispint}{Negative binomial dispersion parameter "intercept" with respect to read counts. See
DESeq2 paper for dispersion model used.}

\item{logkdegsdtrend_slope}{Slope for log10(read count) vs. log(kdeg) replicate variability trend}

\item{logkdegsdtrend_intercept}{Intercept for log10(read count) vs. log(kdeg) replicate variability trend}

\item{logksynsdtrend_slope}{Slope for log10(read count) vs. log(ksyn) replicate variability trend}

\item{logksynsdtrend_intercept}{Intercept for log10(read count) vs. log(ksyn) replicate variability trend}

\item{logkdeg_mean}{Mean of normal distribution from which reference log(kdeg)
linear model parameter is drawn from for each feature if \code{param_details} is not provided.}

\item{logkdeg_sd}{Standard deviation of normal distribution from which reference log(kdeg)
linear model parameter is drawn from for each feature if \code{param_details} is not provided.}

\item{logksyn_mean}{Mean of normal distribution from which reference log(ksyn)
linear model parameter is drawn from for each feature if \code{param_details} is not provided.}

\item{logksyn_sd}{Standard deviation of normal distribution from which reference log(ksyn)
linear model parameter is drawn from for each feature if \code{param_details} is not provided.}

\item{logkdeg_diff_avg}{Mean of normal distribution from which non-reference log(kdeg)
linear model parameters are drawn from for each feature if \code{param_details} is not provided.}

\item{logksyn_diff_avg}{Mean of normal distribution from which reference log(ksyn)
linear model parameter are drawn from for each feature if \code{param_details} is not provided.}

\item{logkdeg_diff_sd}{Standard deviation of normal distribution from which reference log(kdeg)
linear model parameter are drawn from for each feature if \code{param_details} is not provided.}

\item{logksyn_diff_sd}{Standard deviation of normal distribution from which reference log(ksyn)
linear model parameter are drawn from for each feature if \code{param_details} is not provided.}

\item{pdiff_kd}{Proportion of features for which non-reference log(kdeg) linear model parameters
differ from the reference.}

\item{pdiff_ks}{Proportion of features for which non-reference log(ksyn) linear model parameters
differ from the reference.}

\item{pdiff_both}{Proportion of features for which BOTH non-reference log(kdeg) and log(ksyn) linear model parameters
differ from the reference.
ksyns are simulated}

\item{pdo}{Dropout rate; think of this as the probability that a s4U containing
molecule is lost during library preparation and sequencing. If \code{pdo} is 0 (default)
then there is not dropout.}

\item{dynamics_preset}{Which preset model to use for simulation of dynamics.
Therefore, only relevant if \code{mode} == \code{dynamics}. Options are:
\describe{
\item{nuc2cyto}{Simplest model of nuclear and cytoplasmic RNA dynamics: 0 -> N -> C -> 0}
\item{preRNA}{Simplest model of pre-RNA and mature RNA dynamics: 0 -> P -> M -> 0}
\item{preRNAwithPdeg}{Same as preRNA, but now pre-RNA can also degrade.}
\item{nuc2cytowithNdeg}{Same as nuc2cyto, but now nuclear RNA can also degrade.}
\item{subtlseq}{Subcellular TimeLapse-seq model, similar to that described in Ietswaart et al., 2024.
Simplest model discussed there, lacking nuclear degradation: 0 -> CH -> NP -> CY -> PL -> 0, and CY can
also degrade.}
\item{nuc2cytowithpreRNA}{Combination of nuc2cyto and preRNA where preRNA is first synthesized,
then either processed or exported to the cytoplasm. Processing can also occur in the cytoplasm, and
mature nuclear RNA can be exported to the cytoplasm. Only  mature RNA degrades.}
}}

\item{unassigned_name}{String to give to reads not assigned to a given feature.}

\item{dispersion}{Negative binomial \code{size} parameter to use for simulating read counts}

\item{lfn_sd}{Logit(fn) replicate variability.}

\item{treatment_effects}{Data frame describing effects of treatment on each
parameter. Should have five columns: "parameter_index", "treatment_index", "mean", "sd",
and "fraction_affected".
Each row corresponds to the effect the ith (i = treatment_index) treatment has on the
jth (j = parameter_index) kinetic parameter. Effect sizes, on a log-scale, are drawn from
a Normal distribution with mean and standard deviation set by the mean and sd columns,
respectively. The number of non-zero effects is set by "fraction_affected", and is
equal to \code{ceiling(nfeatures * fraction_affected)}. treatment_index of 1 will be ignored
and can either be included or not.}

\item{effect_avg_default}{If \code{ntreatments} > 1, and \code{treatment_effects} is not
provided, this will be the value of \code{mean} for all treatments and parameters imputed
in \code{treatment_effects}.}

\item{effect_sd_default}{If \code{ntreatments} > 1, and \code{treatment_effects} is not
provided, this will be the value of \code{sd} for all treatments and parameters imputed
in \code{treatment_effects}.}

\item{fraction_affected_default}{If \code{ntreatments} > 1, and \code{treatment_effects} is not
provided, this will be the value of \code{fraction_affected} for all treatments and parameters imputed
in \code{treatment_effects}.}

\item{log_means}{Vector of log-Normal logmeans from which the distribution of
feature-specific parameters will be drawn from. Length of vector should be the same
as max(entries in \code{graph}), i.e., the number of parameters in your specified model.
If not provided, will by default be \code{c(1, seq(from = -0.3, to = -2.5, length.out = max(graph) - 1 ))}.
\code{1} for the ksyn parameter (which is always denoted 1 in the preset \code{graph}) is arbitrary.
Remaining parameters will make it so indices order parameters from fastest to slowest process.}

\item{log_sds}{Vector of log-Normal logsds from which the distribution of
feature-specific parameters will be drawn from. If not provided, will be 0.4 for all parameters.}
}
\description{
\code{EZSimulate()} is a user friendly wrapper to \code{SimulateMultiCondition()}. It
sets convenient defaults so as to quickly generate easy to interpret output.
\code{EZSimulate()} has all of the same parameters as \code{SimulateMultiCondition()},
but it also has a number of additional parameters that guide its default behavior
and allow you to simulate multi-condition data without specifying the multiple,
sometimes complex, arguments that you would need to specify in \code{SimulateMultiCondition()}
to get the same behavior. In particular, users only have to set a single parameter,
\code{nfeatures} (number of features to simulate data for), by default. The \code{EZSimulate()}-unique
parameters \code{ntreatments} and \code{nreps} have default values that guide the simulation in the
case where only \code{nfeatures} is specified. In particular, \code{nreps} of \code{ntreatments} different
conditions will be simulated, with the assumed model \code{log(kdeg) ~ treatment} and \code{log(ksyn) ~ 1}.
In other words, Different kdeg values will be simulated for each treatment level, and ksyn
values will not differ across conditions.
}
\examples{

# Simulate standard data
simdata_standard <- EZSimulate(30)

# Simulate dynamical systems data
simdata_ode <- EZSimulate(30,
                          mode = "dynamics",
                          ntreatments = 1,
                          label_time = c(1, 3),
                          dynamics_preset = "nuc2cyto")

}
