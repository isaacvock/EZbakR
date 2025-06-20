% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fractions.R
\name{EstimateMutRates}
\alias{EstimateMutRates}
\alias{EstimateMutRates.EZbakRData}
\alias{EstimateMutRates.EZbakRArrowData}
\title{Estimate mutation rates}
\usage{
EstimateMutRates(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)

\method{EstimateMutRates}{EZbakRData}(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)

\method{EstimateMutRates}{EZbakRArrowData}(
  obj,
  populations = "all",
  pnew_prior_mean = -2.94,
  pnew_prior_sd = 0.3,
  pold_prior_mean = -6.5,
  pold_prior_sd = 0.5,
  pold_est = NULL,
  pold_from_nolabel = FALSE,
  grouping_factors = NULL
)
}
\arguments{
\item{obj}{An \code{EZbakRData} or \code{EZbakRArrowData} object}

\item{populations}{Character vector of the set of mutational populations
that you want to infer the fractions of. For example, say your cB file contains
columns tracking T-to-C and G-to-A}

\item{pnew_prior_mean}{logit-Normal mean for logit(pnew) prior.}

\item{pnew_prior_sd}{logit-Normal sd for logit(pnew) prior.}

\item{pold_prior_mean}{logit-Normal mean for logit(pold) prior.}

\item{pold_prior_sd}{logit-Normal sd for logit(pold) prior.}

\item{pold_est}{Background mutation rate estimates if you have them. Can either be a single
number applied to all samples or a named vector of values, where the names should be sample
names.}

\item{pold_from_nolabel}{Fix background mutation rate estimate to mutation rates seen in -label samples.
By default, a single background rate is used for all samples, inferred from the average mutation rate
across all -label samples. The \code{grouping_factors} argument can be specified to use certain -label samples
to infer background mutation rates for certain sets of +label samples.}

\item{grouping_factors}{If \code{pold_from_nolabel} is TRUE, then \code{grouping_factors} will specify the
sample-detail columns in the metadf that should be used to group -label samples by. Average mutation
rates in each group of -label samples will be used as the background mutation rate estimate in
+label samples with the same values for the relevant metadf columns.}
}
\value{
\code{EZbakRData} object with an added \code{mutrates} slot containing estimated
high and low mutation rates for each mutation type modeled.
}
\description{
Two component mixture models are fit to all data to estimate global high and
low mutation rates for all samples. Estimation of these mutation rates
are regularized through the use of weakly informative priors whose parameters can
be altered using the arguments defined below.
}
\section{Methods (by class)}{
\itemize{
\item \code{EstimateMutRates(EZbakRData)}: Method for class \strong{EZbakRData}
Estimates mutation rates using a fully in memory object.

\item \code{EstimateMutRates(EZbakRArrowData)}: Method for class \strong{EZbakRArrowData}
Estimate mutation rates using a partially on-disk object.

}}
\examples{

# Simulate data to analyze
simdata <- SimulateOneRep(30)

# Create EZbakR input
metadf <- data.frame(sample = "sampleA", tl = 2)
ezbdo <- EZbakRData(simdata$cB, metadf)

# Estimate mutation rates
mutrates <- EstimateMutRates(ezbdo)

}
