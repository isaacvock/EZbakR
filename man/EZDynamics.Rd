% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODEs.R
\name{EZDynamics}
\alias{EZDynamics}
\title{Generalized dynamical systems modeling}
\usage{
EZDynamics(
  obj,
  graph,
  sub_features,
  grouping_features,
  scale_factors = NULL,
  sample_feature = NULL,
  modeled_to_measured = NULL,
  parameter_names = paste0("logk", 1:max(graph)),
  unassigned_name = "__no_feature",
  type = "averages",
  prior_means = rep(-3, times = max(graph)),
  prior_sds = c(3, rep(1, times = max(graph) - 1)),
  avg_params_tokeep = NULL,
  avg_params_todrop = NULL,
  label_time_name = "tl",
  features = NULL,
  populations = NULL,
  fraction_design = NULL,
  parameter = NULL,
  repeatID = NULL,
  exactMatch = TRUE,
  feature_lengths = NULL,
  use_coverage = TRUE,
  normalization_repeatID = NULL,
  normalization_exactMatch = TRUE,
  species_to_sf = NULL,
  overwrite = TRUE
)
}
\arguments{
\item{obj}{Currently must be an EZbakRData object on which \code{AverageAndRegularize}
has been run. In the future, will also support (in case where all species are
assayed in every sample) providing output of just \code{EstimateFractions()} as input,
acting as a generalization of \code{EstimateKinetics()} in that case.}

\item{graph}{An NxN adjacency matrix, where N represents the number of species
being modeled. One of these species must be called "0" and represent the "no RNA"
species. This is the species from which some species are synthesized (e.g., 0 -> P,
means premature RNA is synthesized from no RNA), and the species to which
some species are degraded (e.g., M -> 0 means mature RNA is converted to "no RNA"
via degradation). The rows and columns of this matrix must be the names of all
modled species, and rownames(graph) == colnames(graph). Entry i,j of the matrix
is either 0 if species i cannot be converted into species j under your model,
and an integer from 1:npars (where npars = total number of parameters to be
estimated) if it can.

For example, the model 0 -> P -> M -> 0 would have the \code{graph}:
\verb{matrix(c(0, 1, 0, 0, 0, 2, 3, 0, 0), nrow = 3, ncol = 3, byrow = TRUE}.}

\item{sub_features}{Which feature columns distinguish between the different
measured species? Note, the measured species need not have the same name,
and may not be directly equivalent to, the modeled species. The relationship
between the modeled species in \code{graph} and \code{sub_features} needs to be specified
in \code{modeled_to_measured} if the names are not equivalent though.}

\item{grouping_features}{Which features are the overarching feature assignments
by which \code{sub_features} should be grouped? This will usually be the feature
columns corresponding to full-gene assignments, as well as any higher order
assignments (e.g., chromosome). A \code{sub_feature} can be included in \code{grouping_features}
if it never has the value of \code{unassigned_name} ("__no_feature" by default). Only
one \code{sub_feature} should ever fulfill this criterion though.}

\item{scale_factors}{Data frame mapping samples to factors by which to multiply
read counts so as ensure proper normalization between different RNA populations.
Only relevant if you are modeling relationships between distinct RNA populations,
for example RNA from nuclear and cytoplasmic fractions. Will eventually be inferred
automatically.}

\item{sample_feature}{If different samples involve assaying different species,
this must be the name of the metadf column that distinguishes the different
classes of samples. For example, if analyzing a subcellular fractionation
dataset, you likely included a column in your metadf called "compartment".
This would then be your \code{sample_feature}, assuming you ran \code{AverageAndRegularize()},
with a \code{mean_formula} that included compartment as a term.}

\item{modeled_to_measured}{If \code{sub_features} is not identical to the non "0"
species names in \code{graph}, then you must specify the relationship between
\code{sub_features}, \code{sample_feature} (if specified), and the species in \code{graph}.
This is done through a list of formulas, whose left side is a \code{sub_feature}
and whose right side is a function of species in \code{graph}. If \code{sample_feature}
is not specified, then \code{modeled_to_measured} should be a nameless list of
formulae. If \code{sample_feature} is specified, then \code{modeled_to_measured} must
be a named list of formulas where the names correspond to unique values of
\code{sample_feature}. In this latter case, the elements, should be the mapping
of measured features (\code{sub_features}) to modeled species (names in \code{graph}
that aren't "0").

For example, if your model is 0 -> P -> M -> 0, where P stands for premature RNA
and M stands for mature RNA, and you have a column called
GF that corresponds to assignment anywhere in the gene, and XF that corresponds
to assignment of exclusively exon mapping reads, then your \code{modeled_to_measured}
should be \verb{list(GF ~ P + M, XF ~ M).}.

As another example, if your model is 0 -> N -> C -> 0, where N stands for
"nuclear RNA" and C stands for "cytoplasmic RNA", and your \code{sample_feature}
takes on values of "nuclear", "cytoplasmic", and "total", and you have a single
\code{sub_feature} called XF, then your \code{modeled_to_measured} shouild be
\code{list(nuclear = GF ~ N, cytoplasmic = GF ~ C, total = GF ~ C + N)}. This is
interpreted as meaning in groups of samples for which \code{sample_feature} == "nuclear",
reads assigned to a GF are from the N species (nuclear RNA). When
\code{sample_feature} == "cytoplasmic", reads assigned to GF correspond to the C
species (cytoplasmic RNA). When \code{sample_feature} == "total", reads assigned
to GF correspond to a combination of N and C (nuclear and cytoplasmic RNA).}

\item{parameter_names}{Vector of names you would like to give to the estimated
parameters. ith element should correspond to name of parameter given the ID
i in \code{graph}. By default, this is just ki, where i is this numerical index.}

\item{unassigned_name}{What value will a \code{sub_feature} column have if a read
was not assigned to said feature? "__no_feature" by default.}

\item{type}{What type of table would you like to use? Currently only supports
"averages", but will support "fractions" in the near future.}

\item{prior_means}{Mean of log-Normal prior for kinetic parameters. Should
be vector where ith value is mean for ith parameter, i = index in \code{graph}}

\item{prior_sds}{Std. dev. of log-Normal prior for kinetic parameter.
Should be vector where ith value is mean for ith parameter, i = index in \code{graph}.}

\item{avg_params_tokeep}{Names of parameters in averages table that you would
like to keep. Other parameters will be discarded. Don't include the prefixes "mean_",
"sd_", or "coverage_"; these will be imputed automatically. In other words,
this should be the base parameter names.}

\item{avg_params_todrop}{Names of parameters in averages table that you would
like to drop. Other parameters will be kept. Don't include the prefixes "mean_",
"sd_", or "coverage_"; these will be imputed automatically. In other words,
this should be the base parameter names.}

\item{label_time_name}{Name of relevant label time column that will be found
in the parameter names. Defaults to the standard "tl".}

\item{features}{Character vector of the set of features you want to stratify
reads by and estimate proportions of each RNA population. The default of "all"
will use all feature columns in the \code{obj}'s cB.}

\item{populations}{Mutational populations that were analyzed to generate the
fractions table to use. For example, this would be "TC" for a standard
s4U-based nucleotide recoding experiment.}

\item{fraction_design}{"Design matrix" specifying which RNA populations exist
in your samples. By default, this will be created automatically and will assume
that all combinations of the \code{mutrate_populations} you have requested to analyze are
present in your data. If this is not the case for your data, then you will have
to create one manually. See docs for \code{EstimateFractions} (run ?EstimateFractions()) for more details.}

\item{parameter}{Parameter to average across replicates of a given condition.
Has to be "logit_fraction_high\if{html}{\out{<muttype>}}", where \if{html}{\out{<muttype>}} is the type of
mutation modeled in \code{EstimateFractions()} (e.g, TC) in this case.}

\item{repeatID}{If multiple \code{fractions} tables exist with the same metadata,
then this is the numerical index by which they are distinguished.}

\item{exactMatch}{If TRUE, then \code{features} and \code{populations} have to exactly match
those for a given fractions table for that table to be used. Means that you can't
specify a subset of features or populations by default, since this is TRUE
by default.}

\item{feature_lengths}{Table of effective lengths for each feature combination in your
data. For example, if your analysis includes features named GF and XF, this
should be a data frame with columns GF, XF, and length.}

\item{use_coverage}{If TRUE, normalized read counts will be used to inform
kinetic parameter estimates. If FALSE, only fraction news will be used, which
will leave some parameters (e.g., synthesis rate) unidentifiable, though has
the advantage of avoiding the potential biases induced by normalization problems.}

\item{normalization_repeatID}{For extracting the \code{fractions} table needed for
normalization of multi-sample data. If multiple \code{fractions} tables exist with the same metadata,
then this is the numerical index by which they are distinguished.}

\item{normalization_exactMatch}{For extracting the \code{fractions} table needed for
normalization of multi-sample data. If TRUE, then \code{features} and \code{populations}
have to exactly match those for a given fractions table for that table to be used.
Means that you can't specify a subset of features or populations by default, since this is TRUE
by default.}

\item{species_to_sf}{List mapping individual RNA species in \code{graph} to
different sample_feature values (sf). This is relevant if you are modeling both
pre- and mature RNA dynamics in subcellular fractionation data. EZbakR can usually
automatically infer this, but if not, then you can manually specify this mapping.
Should be a list with one element per unique sample_feature type. Each element should
be a vector of the RNA species in \code{graph} that belong to that sample_feature type. For
example, if you have whole cell, cytoplasmic, and nuclear fraction data, this should have
one element called "cytoplasmic" and one element called "nuclear". The whole cell
sample_feature is a sum of cytoplasmic and nuclear and thus does not apply. The "cytoplasmic"
element of this list should be the set of RNA species that are present in the cytoplasm,
e.g. cytoplasmic pre-RNA (maybe referred to in \code{graph} as CP) and cytoplasmic mature
RNA (maybe referred to in \code{graph} as CM).}

\item{overwrite}{If TRUE and a fractions estimate output already exists that
would possess the same metadata (features analyzed, populations analyzed,
and fraction_design), then it will get overwritten with the new output. Else,
it will be saved as a separate output with the same name + "_#" where "#" is a
numerical ID to distinguish the similar outputs.}
}
\value{
\code{EZbakRData} object with an additional "dynamics" table.
}
\description{
\code{EZDynamics()} estimates parameters of a user-specified dynamical systems
model. The dynamical system model is specified through an adjacency matrix,
which is an NxN matrix described below (see \code{graph} documentation). Modeling
can either be done for species all assayed in each sample, or species that
are assayed across a set of independent samples (e.g., subcellular fractionation
involves assaying different species in different samples).
}
\details{
When running \code{AverageAndRegularize()} to produce input for \code{EZDynamics()}, you
must set \code{parameter} to "logit_fraction_high\if{html}{\out{<muttype>}}" (\if{html}{\out{<muttype>}} = type of mutation
modeled by \code{EstimateFractions()}, e.g., TC). If you have multiple distinct label times,
you must also include the label time (\code{tl} of your \code{metadf})
in your regression formula. \code{EZDynamics()} models the logit(fraction high <muttype),
and this will depend on the label time (longer label time = higher fraction), which
is why these two conditions must be met. If you only have a single label time though,
\code{EZDynamics} will be able to impute this one value for all samples from your \code{metadf}.
You can also include additional interaction terms in your \code{AverageAndRegularize()}
model for different experimental conditions in which experiments were conducted,
so that inferred kinetic parameters can be compared across these conditions. Currently,
more complex modeling beyond simple stratification of samples by one or more condition
is not possible with \code{EZDynamics()}.

For normalization purposes, especially if analyzing pre-mRNA processing dynamics,
you will need to provide \code{AverageAndRegularize()} with a table of feature lengths
via the \code{feature_lengths} parameter. This will be used in all cases to length
normalize read counts. Even in the case when you are just modeling mature mRNA
dynamics, this is technically necessary for accurate estimation of scale factors.

The first step of \code{EZDynamics()} is attempted inference of normalization scale
factors for read counts. If you have scale factors you calculated yourself,
e.g. via specialized spike-in protocols, you can provide these via the \code{scale_factors}
parameter. If not, \code{EZDynamics()} will try to infer these from the fraction labeled's
in each \code{sample_feature} (e.g., in different subcellular compartments). This relies on
having some \code{sample_feature}'s that are a combination of other \code{sample_feature}'s. For
example, if analyzing subcellular fractionation data, you may have 1) total RNA; 2)
cytoplasmic RNA; and 3) nuclear RNA. Total RNA = cytoplasmic + nuclear RNA, and thus
the fraction of reads from labeled RNA is a function of the total cytoplasmic and
nuclear fraction labeled's, as well as the relative molecular abundances of cytoplasmic
and nuclear RNA. The latter is precisely the scale factors we need to estimate.
If you do not have sufficient combinations of data to perform this scale factor estimation,
\code{EZDynamics()} will only use the fraction labeled's for modeling kinetic parameters.
It can then perform post-hoc normalization to estimate a single synthesis rate constant,
using the downstream rate constants to infer the unknown normalization scale factor necessary
to combine kinetic parameter estimates and read counts to infer this rate constant.

For estimating kinetic parameters, \code{EZDynamics()} infers the solution of the linear
system of ODEs specified in your \code{graph} matrix input. This is done by representing
the system of equations as a matrix, and deriving the general solution of the levels
of each modeled species of RNA from the eigenvalues and eigenvectors of this matrix.
While this makes \code{EZDynamics()} orders of magnitude more efficient than if it had to
numerically infer the solution for each round of optimization, needing to compute
eigenvalues and eigenvectors in each optimization iteration is still non-trivial,
meaning that \code{EZDynamics()} may take anywhere from 10s of minutes to a couple hours
to run, depending on how complex your model is and how many distinct set of samples
and experimental conditions you have.
}
\examples{

##### MODELING CYTOPLASMIC TO NUCLEAR FLOW

### Simulate data and get replicate average fractions estimates
simdata <- EZSimulate(
  nfeatures = 30,
  ntreatments = 1,
  mode = "dynamics",
  label_time = c(1, 3),
  dynamics_preset = "nuc2cyto"
)

ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
ezbdo <- EstimateFractions(ezbdo)
ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~tl:compartment - 1,
                              type = "fractions",
                              parameter = "logit_fraction_highTC")

### ODE model: the graph and measured species
graph <- matrix(c(0, 1, 0,
                  0, 0, 2,
                  3, 0, 0),
                nrow = 3,
                ncol = 3,
                byrow = TRUE)
colnames(graph) <- c("0", "N", "C")
rownames(graph) <- colnames(graph)

modeled_to_measured <- list(
  nuclear = list(GF ~ N),
  cytoplasm = list(GF ~ C),
  total = list(GF ~ C + N) # total RNA is a combination of C and N
)

### Fit model
ezbdo <- EZDynamics(ezbdo,
                    graph = graph,
                    sub_features = "GF",
                    grouping_features = "GF",
                    sample_feature = "compartment",
                    modeled_to_measured = ode_models$nuc2cyto$formulas)


}
