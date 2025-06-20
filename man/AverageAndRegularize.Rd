% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Comparisons.R
\name{AverageAndRegularize}
\alias{AverageAndRegularize}
\title{Average parameter estimates across replicates, and regularize variance estimates}
\usage{
AverageAndRegularize(
  obj,
  features = NULL,
  parameter = "log_kdeg",
  type = "kinetics",
  kstrat = NULL,
  populations = NULL,
  fraction_design = NULL,
  exactMatch = TRUE,
  repeatID = NULL,
  formula_mean = NULL,
  sd_grouping_factors = NULL,
  include_all_parameters = TRUE,
  sd_reg_factor = 10,
  error_if_singular = TRUE,
  min_reads = 10,
  convert_tl_to_factor = TRUE,
  regress_se_with_abs = TRUE,
  force_lm = FALSE,
  force_optim = force_lm,
  conservative = FALSE,
  character_limit = 20,
  feature_lengths = NULL,
  feature_sample_counts = NULL,
  scale_factor_df = NULL,
  overwrite = TRUE
)
}
\arguments{
\item{obj}{An \code{EZbakRFractions} or \code{EZbakRKinetics} object, which is an \code{EZbakRData} object on
which \code{EstimateFractions()} or \code{EstimateKinetics()} has been run.}

\item{features}{Character vector of the set of features you want to stratify
reads by and estimate proportions of each RNA population. The default of "all"
will use all feature columns in the \code{obj}'s cB.}

\item{parameter}{Parameter to average across replicates of a given condition.}

\item{type}{What type of table is the parameter found in? Default is "kinetics",
but can also set to "fractions".}

\item{kstrat}{If \code{type == "kinetics"}, then \code{kstrat} specifies the kinetic parameter
inference strategy.}

\item{populations}{Character vector of the set of mutational populations
that you want to infer the fractions of. Only relevant if type == "fractions".}

\item{fraction_design}{"Design matrix" specifying which RNA populations exist
in your samples. Only relevant if type == "fractions".}

\item{exactMatch}{If TRUE, then \code{features} and \code{populations} have to exactly match
those for a given fractions table for that table to be used. Means that you can't
specify a subset of features or populations by default, since this is TRUE
by default.}

\item{repeatID}{If multiple \code{kinetics} or \code{fractions} tables exist with the same metadata,
then this is the numerical index by which they are distinguished.}

\item{formula_mean}{An R formula object specifying how the \code{parameter} of interest
depends on the sample characteristics specified in \code{obj}'s metadf. The most common formula
will be \code{~ treatment} or \code{~ treatment:duration}, where \code{treatment} and \code{duration} would
be replaced with whatever you called the relevant sample characteristics in your metadf.
\code{~ treatment} means that an average value of \code{parameter} should be estimated for each
set of samples with the same value for \code{treatment} in the metadf. \code{~ treatment:duration} specifies
that an average value of \code{parameter} should be estimated for each set of samples with the same
combination of \code{treatment} and \code{duration} values in the metadf. An example of the latter
case is a situation where you have two or more treatments (e.g., drug treated and untreated control)
which were applied for different durations of time (e.g., 4 and 8 hours).

NOTE: EZbakR automatically removes any intercept terms from the model. That way,
there is no ambiguity about what parameter is defined as the reference.}

\item{sd_grouping_factors}{What metadf columns should data be grouped by when estimating
standard deviations across replicates? If this is NULL, then EZbakR will check to see
if the \code{formula_mean} specifies a formula that cleanly stratifies samples into disjoint
groups. For example, the formula \code{~ treatment} will assign each sample to a single factor
(its value for the metadf's \code{treatment} column). In this case, standard deviations can
be calculated for sets of replicates in each \code{treatment} group. If such a stratification
does not exist, a single standard deviation will be estimated for each feature (i.e.,
homoskedasticity will be assumed as in standard linear modeling).}

\item{include_all_parameters}{If TRUE, an additional table will be saved with the prefix \code{fullfit_},
which includes all of the parameters estimated throughout the course of linear modeling and
regularization. This can be nice for visualizing the regularized mean-variance trend.}

\item{sd_reg_factor}{Determines how strongly variance estimates are shrunk towards trend.
Higher numbers lead to more regularization. Eventually, this will be replaced with estimation
of how much variance there seems to be in the population of variances.}

\item{error_if_singular}{If TRUE, linear model will throw an error if parameters
cannot be uniquely identified. This is most often caused by parameters that cannot
be estimated from the data, e.g., due to limited replicate numbers or correlated
sample characteristics (i.e., all treatment As also correspond to batch As, and
all treatment Bs correspond to batch Bs).}

\item{min_reads}{Minimum number of reads in all samples for a feature to be kept.}

\item{convert_tl_to_factor}{If a label time variable is included in the \code{formula_mean},
convert its values to factors so as to avoid performing continuous regression on label
times. Defaults to TRUE as including label time in the regression is often meant to
stratify samples by their label time if, for example, you are averaging logit(fractions).}

\item{regress_se_with_abs}{If TRUE, and if \code{type == "fractions"}, then standard error
will be regressed against logit fraction rather than magnitude of logit fraction.
Makes sense to set this to FALSE if analyzing certain site-specific mutational probing
methods when high mutation content things are likely low variance SNPs.}

\item{force_lm}{Certain formula lend them selves to efficient approximations of the
full call to \code{lm()}. Namely, formulas that stratify samples into disjoint groups where
a single parameter of the model is effectively estimated from each group can be tackled
via simple averaging of data from each from group. If you would like to force EZbakR
to fit the fully rigorous linear model though, set \code{force_lm} to \code{TRUE}.}

\item{force_optim}{Old parameter that is now passed the value \code{force_lm} and will be
deprecated in later releases}

\item{conservative}{If TRUE, conservative variance regularation will be performed.
In this case, variances below the trend will be regularized up to the trend, and
variances above the trend will be left unregularized. This avoids underestimation
of variances.}

\item{character_limit}{Limit on the number of characters of the name given to the
output table. Will attempt to concatenate the parameter name with the names of all
of the features. If this is too long, only the parameter name will be used.}

\item{feature_lengths}{Table of effective lengths for each feature combination in your
data. For example, if your analysis includes features named GF and XF, this
should be a data frame with columns GF, XF, and length.}

\item{feature_sample_counts}{Data frame with columns \if{html}{\out{<feature names>}} and \code{nsamps},
where \if{html}{\out{<feature names>}} are all of the feature columns in the input to
\code{AverageAndRegularize()}, and \code{nsamps} is the number of samples that samples
from that feature combination needs to have over the read count threshold.}

\item{scale_factor_df}{Data frame with columns "sample" and a second column of
whatever name you please. The second column should denote scale factors by which
read counts in that sample should be multiplied by in order to normalize these
read counts.}

\item{overwrite}{If TRUE, identical, existing output will be overwritten.}
}
\value{
\code{EZbakRData} object with an additional "averages" table, as well as a
fullfit table under the same list heading, which includes extra information about
the priors used for regularization purposes.
}
\description{
\code{AverageAndRegularize} fits a generalized linear model to your data to effectively
average parameter estimates across replicates and get overall uncertainty estimates
for those parameters. The linear model to which your data is fit is specified via
an R formula object supplied to the \code{formula_mean} parameter. Uncertainty estimates
are regularized via a hierarchical modeling strategy originally introduced with
bakR, though slightly improved upon since then.
}
\details{
The EZbakR website has an extensive vignette walking through various use cases
and model types that you can fit with \code{AverageAndRegularize()}: \href{https://isaacvock.github.io/EZbakR/articles/Linear-modeling.html}{vignette link}.
EZbakR improves upon bakR by balancing extra conservativeness in several steps
with a more highly powered statistical testing scheme in its \code{CompareParameters()}
function. In particular, the following changes to the variance regularization
scheme were made:
\itemize{
\item Sample-specific parameter uncertainties are used to generate conservative estimates
of feature-specific replicate variabilties. In addition, a small floor is set to ensure
that replicate variance estimates are never below a certain level, for the same reason.
\item Condition-wide replicate variabilities are regressed against both read coverage and
either a) |logit(estimate)| when modeling average fraction labeled. This captures
the fact thta estimates are best around a logit(fraction labeled) of 0 and get
worse for more extreme fraction labeled's.; b) log(kdeg) when modeling log degradation
rate constants. At first, I considered a strategy similar to the fraction labeled
modeling, but found that agreement between a fully rigorous MCMC sampling approach
and EZbakR was significantly improved by just regressing hee value of the log kientic parameter,
likely due to the non-linear transformation of fraction labeled to log(kdeg);
and c) only coverage in all other cases.
\item Features with replicate variabilities below the inferred trend have their replicate
variabilites set equal to that predicted by the trend. This helps limit underestimation
of parameter variance. Features with above-trend replicate variabilties have their
replicate variabilities regularized with a Normal prior Normal likelihood Bayesian
model, as in bakR (so the log(variance) is the inferred mean of this distribution, and
the known variance is inferred from the amount of variance about the linear dataset-wide
trend).
}
All of this allows \code{CompareParameters()} to use a less conservative statistical test
when calculating p-values, while still controlling false discovery rates.
}
\examples{

# Simulate data to analyze
simdata <- EZSimulate(30)

# Create EZbakR input
ezbdo <- EZbakRData(simdata$cB, simdata$metadf)

# Estimate Fractions
ezbdo <- EstimateFractions(ezbdo)

# Estimate Kinetics
ezbdo <- EstimateKinetics(ezbdo)

# Average estimates across replicate
ezbdo <- AverageAndRegularize(ezbdo)

}
