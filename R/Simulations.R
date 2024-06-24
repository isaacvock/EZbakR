#' Simulate a single replicate of NR-seq data
#'
#' In `SimulateOneRep`, users have the option to either provide vectors of feature-specific
#' read counts, fraction news, kdegs, and ksyns for the simulation, or to have those drawn
#' from relevant distributions whose properties can be tuned by the various optional
#' parameters of `SimulateOneRep`. The number of mutable nucleotides (nT) in
#' a read is drawn from a binomial distribution with `readlength` trials and a probability
#' of "success" equal to `Ucont`. A read's status as new or old is drawn from a Bernoulli
#' distribution with probability of "success" equal to the feature's fraction new. If a read
#' is new, the number of mutations in the read is drawn from a binomial distribution with
#' probability of mutation equal to pnew. If a read is old, the number of mutations is instead
#' drawn from a binomial distribution with probability of mutation equal to pold.
#'
#' @param nfeatures Number of "features" (e.g., genes) to simulate data for
#' @param read_vect Vector of length = `nfeatures`; specifies the number of reads
#' to be simulated for each feature. If this is not provided, the number of reads
#' simulated is equal to `round(seqdepth * (ksyn_i/kdeg_i)/sum(ksyn/kdeg))`. In other words,
#' the normalized steady-state abundance of a feature is multiplied by the total number
#' of reads to be simulated and rounded to the nearest integer.
#' @param label_time Length of s^{4}U feed to simulate.
#' @param sample_name Character vector to assign to `sample` column of output simulated
#' data table (the cB table).
#' @param feature_prefix Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows up in the
#' `feature` column of the output simulated data table.
#' @param fn_vect Vector of length = `nfeatures`; specifies the fraction new to use for each
#' feature's simulation. If this is not provided and `kdeg_vect` is, then `fn_vect = 1 - exp(-kdeg_vect*label_time)`.
#' If both `fn_vect` and `kdeg_vect` are not provided, then kdegs are simulated from a joint distribution as
#' described below and converted to a `fn_vect` as when `kdeg_vect` is user-provided.
#' @param kdeg_vect Vector of length = `nfeatures`; specifies the degradation rate constant to use for each
#' feature's simulation. If this is not provided and `fn_vect` is, then `kdeg_vect = -log(1 - fn_vect)/label_time`.
#' If both `kdeg_vect` and `fn_vect` are not provided, each feature's `kdeg_vect` value is drawn from a log-normal distrubition
#' with meanlog = `logkdeg_mean` and sdlog = `logkdeg_sd`. `kdeg_vect` is actually only simulated in the case
#' where `read_vect` is also not provided, as it will be used to simulate read counts as described above.
#' @param ksyn_vect Vector of length = `nfeatures`; specifies the synthesis rate constant to use for each
#' feature's simulation. If this is not provided, and `read_vect` is also not provided, then each
#' feature's `ksyn_vect` value is drawn from a log-normal distribution with meanlog = `logksyn_mean` and
#' sdlog = `logksyn_sd`. ksyn's do not need to be simulated if `read_vect` is provided, as they only
#' influence read counts.
#' @param pnew Probability that a T is mutated to a C if a read is new.
#' @param pold Probability that a T is mutated to a C if a read is old.
#' @param logkdeg_mean If necessary, meanlog of a log-normal distribution from which
#' kdegs are simulated
#' @param logkdeg_sd If necessary, sdlog of a log-normal distribution from which
#' kdegs are simulated
#' @param logksyn_mean If necessary, meanlog of a log-normal distribution from which
#' ksyns are simulated
#' @param logksyn_sd If necessary, sdlog of a log-normal distribution from which
#' ksyns are simulated
#' @param seqdepth Only relevant if `read_vect` is not provided; in that case, this is
#' the total number of reads to simulate.
#' @param readlength Length of simulated reads. In this simple simulation, all reads
#' are simulated as being exactly this length.
#' @param Ucont Probability that a nucleotide in a simulated read is a U.
#' @param feature_pnew Boolean; if TRUE, simulate a different pnew for each feature
#' @param pnew_kdeg_corr Boolean; only relevant if `feature_pnew` is TRUE. If so, then
#' setting `pnew_kdeg_corr` to TRUE will ensure that higher kdeg transcripts have a higher
#' pnew.
#' @param logit_pnew_mean If `feature_pnew` is TRUE, then the logit(pnew) for each feature
#' will be drawn from a normal distribution with this mean.
#' @param logit_pnew_sd If `feature_pnew` is TRUE, then the logit(pnew) for each feature
#' will be drawn from a normal distribution with this standard deviation.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
SimulateOneRep <- function(nfeatures, read_vect = NULL, label_time = 2,
                           sample_name = "sampleA",
                           feature_prefix = "Gene",
                           fn_vect = NULL, kdeg_vect = NULL, ksyn_vect = NULL,
                           pnew = 0.05, pold = 0.002,
                           logkdeg_mean = -1.9, logkdeg_sd = 0.7,
                           logksyn_mean = 2.3, logksyn_sd = 0.7,
                           seqdepth = 10000000, readlength = 200,
                           Ucont = 0.25, feature_pnew = FALSE,
                           pnew_kdeg_corr = FALSE,
                           logit_pnew_mean = -2.5, logit_pnew_sd = 0.1){

  ### Hack to deal with devtools::check() NOTEs
  feature <- TC <- nT <- NULL

  `.` <- list


  ### Check validity of input

  args <- c(as.list(environment()))

  check_SimulateOneRep_input(args)


  ### Simulate feature-specific pnew and pold as necessary

  if(feature_pnew){

    pnew <- inv_logit(stats::rnorm(nfeatures,
                                   logit_pnew_mean,
                                   logit_pnew_sd))

  }



  ### Simulate kinetic parameters as necesary

  # kdeg and fraction new

  if(is.null(fn_vect)){

    if(is.null(kdeg_vect)){

      kdeg_vect <- stats::rlnorm(nfeatures,
                          logkdeg_mean,
                          logkdeg_sd)


    }

    fn_vect <- 1 - exp(-kdeg_vect*label_time)


  }

  # read counts

  if(is.null(read_vect)){

    if(is.null(ksyn_vect)){

      ksyn_vect <- stats::rlnorm(nfeatures,
                          logksyn_mean,
                          logksyn_sd)

    }

    if(is.null(kdeg_vect)){

      kdeg_vect <- -log(1 - fn_vect)/label_time

    }

    read_vect <- round(((ksyn_vect/kdeg_vect)/sum(ksyn_vect/kdeg_vect))*seqdepth)

  }

  if(length(read_vect) == 1 & nfeatures > 1){

    read_vect <- rep(read_vect, times = nfeatures)

  }


  ### Simulate mutational data

  totreads <- sum(read_vect)

  read_status <- stats::rbinom(n = totreads,
                        size = 1,
                        prob = rep(fn_vect, times = read_vect))

  nT_count <- stats::rbinom(n = totreads,
                     size = readlength,
                     prob = Ucont)

  if(feature_pnew){

    if(pnew_kdeg_corr){

      pnew <- pnew[order(pnew)]

      TC_count <- stats::rbinom(n = totreads,
                                size = nT_count,
                                prob = read_status*rep(pnew[rank(kdeg_vect)], times = read_vect) + (1 - read_status)*pold)


    }else{

      TC_count <- stats::rbinom(n = totreads,
                                size = nT_count,
                                prob = read_status*rep(pnew, times = read_vect) + (1 - read_status)*pold)


    }



  }else{

    TC_count <- stats::rbinom(n = totreads,
                              size = nT_count,
                              prob = read_status*pnew + (1 - read_status)*pold)

  }


  cB <- data.table::data.table(
    sample = sample_name,
    feature = rep(paste0(feature_prefix, 1:nfeatures),
                  times = read_vect),
    TC = TC_count,
    nT = nT_count
  )[,.(n = .N), by = .(sample, feature, TC, nT)]


  ### Save ground truth

  if(feature_pnew){

    if(pnew_kdeg_corr){

      truth <- data.table::data.table(sample = sample_name,
                                      feature = paste0(feature_prefix, 1:nfeatures),
                                      true_fraction_highTC = fn_vect,
                                      true_kdeg = kdeg_vect,
                                      true_ksyn = ksyn_vect,
                                      true_pnew = pnew[rank(kdeg_vect)])

    }else{

      truth <- data.table::data.table(sample = sample_name,
                                      feature = paste0(feature_prefix, 1:nfeatures),
                                      true_fraction_highTC = fn_vect,
                                      true_kdeg = kdeg_vect,
                                      true_ksyn = ksyn_vect,
                                      true_pnew = pnew)

    }


  }else{

    truth <- data.table::data.table(sample = sample_name,
                                    feature = paste0(feature_prefix, 1:nfeatures),
                                    true_fraction_highTC = fn_vect,
                                    true_kdeg = kdeg_vect,
                                    true_ksyn = ksyn_vect)

  }




  return(list(cB = cB,
              ground_truth = truth))



}


#' Simulate NR-seq data for multiple replicates of multiple biological conditions
#'
#' `EZSimulate()` is a user friendly wrapper to `SimulateMultiCondition()`. It
#' sets convenient defaults so as to quickly generate easy to interpret output.
#' `EZSimulate()` has all of the same parameters as `SimulateMultiCondition()`,
#' but it also has a number of additional parameters that guide its default behavior
#' and allow you to simulate multi-condition data without specifying the multiple,
#' sometimes complex, arguments that you would need to specify in `SimulateMultiCondition()`
#' to get the same behavior. In particular, users only have to set a single parameter,
#' `nfeatures` (number of features to simulate data for), by default. The `EZSimulate()`-unique
#' parameters `ntreatments` and `nreps` have default values that guide the simulation in the
#' case where only `nfeatures` is specified. In particular, `nreps` of `ntreatments` different
#' conditions will be simulated, with the assumed model `log(kdeg) ~ treatment` and `log(ksyn) ~ 1`.
#' In other words, Different kdeg values will be simulated for each treatment level, and ksyn
#' values will not differ across conditions.
#'
#' @param nfeatures Number of "features" (e.g., genes) for which to simulated data.
#' @param ntreatments Number of distinct treatments to simulate. This parameter is
#' only relevant if `metadf` is not provided.
#' @param nreps Number of replicates of each treatment to simulate. This parameter is
#' only relevant if `metadf` is not provided
#' @param metadf A data frame with the following columns:
#' \itemize{
#'  \item sample: Names given to samples to simulate.
#'  \item <details>: Any number of columns with any names (not taken by other metadf columns)
#'  storing factors by which the samples can be stratified. These can be referenced
#'  in `mean_formula`, described below.
#' }
#' These parameters (described more below) can also be included in metadf to specify sample-specific simulation
#' parameter:
#' \itemize{
#'  \item seqdepth
#'  \item label_time
#'  \item pnew
#'  \item pold
#'  \item readlength
#'  \item Ucont
#' }
#' @param mean_formula A formula object that specifies the linear model used to
#' relate the factors in the <details> columns of `metadf` to average log(kdegs) and
#' log(ksyns) in each sample.
#' @param param_details A data frame with one row for each column of the design matrix
#' obtained from `model.matrix(mean_formula, metadf)` that describes how to simulate
#' the linear model parameters. The columns of this data frame are:
#' \itemize{
#'  \item param: Name of linear model parameter as it appears in the column names of the
#'  design matrix from `model.matrix(mean_formula, metadf)`.
#'  \item reference: Boolean; TRUE if you want to treat that parameter as a "reference". This
#'  means that all other parameter values that aren't global parameters are set equal to this
#'  unless otherwise determined (see `pdiff_*` parameters for how it is determined if a parameter
#'  will differ from the reference).
#'  \item global: Boolean; TRUE if you want to treat that parameter as a global parameter. This means
#'  that a single value is used for all features.
#'  \item logkdeg_mean: If parameter is the reference, then its value for the log(kdeg) linear model
#'  will be drawn from a normal distribution with this mean. If it is a global parameter, then this
#'  value will be used. If it is neither of these, then its value in the log(kdeg) linear model will
#'  either be the reference (if there is no difference between this condition's value and the reference)
#'  or the reference's value + a normally distributed random variable centered on this value.
#'  \item logkdeg_sd: sd used for draws from normal distribution as described for `logkdeg_mean`.
#'  \item logksyn_mean: Same as `logkdeg_mean` but for log(ksyn) linear model.
#'  \item logksyn_sd: Same as `logkdeg_sd` but for log(kdeg) linear model.
#'  \item pdiff_ks: Proportion of features whose value of this parameter in the log(ksyn) linear model
#'  will differ from the reference's. Should be a number between 0 and 1, inclusive. For example, if
#'  `pdiff_ks` is 0.1, then for 10% of features, this parameter will equal the reference parameter +
#'  a normally distributed random variable with mean `logksyn_mean` and sd `logksyn_sd`. For the other
#'  90% of features, this parameter will equal the reference.
#'  \item pdiff_kd: Same as `pdiff_ks` but for log(kdeg) linear model.
#'  \item pdiff_both: Proportion of features whose value for this parameter in BOTH the
#'  log(kdeg) and log(ksyn) linear models will differ from the reference. Value must be
#'  between 0 and min(c(pdiff_kd, pdiff_ks)) in that row.
#' }
#' If param_details is not specified by the user, the first column of the design matrix
#' is assumed to represent the reference parameter, all parameters are assumed to be
#' non-global, logkdeg_mean and logksyn_mean are set to the equivalently named parameter values
#' described below for the reference and `logkdeg_diff_avg` and `logksyn_diff_avg` for all other parameters,
#' logkdeg_sd and logksyn_sd are set to the equivalently named parameter values
#' described below for the reference and `logkdeg_diff_sd` and `logksyn_diff_sd` for all other parameters,
#' and pdiff_kd, pdiff_ks, and pdiff_both are all set to the equivalently named parameter values.
#' @param seqdepth Only relevant if `read_vect` is not provided; in that case, this is
#' the total number of reads to simulate.
#' @param label_time Length of s^{4}U feed to simulate.
#' @param pnew Probability that a T is mutated to a C if a read is new.
#' @param pold Probability that a T is mutated to a C if a read is old.
#' @param readlength Length of simulated reads. In this simple simulation, all reads
#' are simulated as being exactly this length.
#' @param Ucont Probability that a nucleotide in a simulated read is a U.
#' @param logkdeg_mean Mean of normal distribution from which reference log(kdeg)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param feature_prefix Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows up in the
#' `feature` column of the output simulated data table.
#' @param dispslope Negative binomial dispersion parameter "slope" with respect to read counts. See
#' DESeq2 paper for dispersion model used.
#' @param dispint Negative binomial dispersion parameter "intercept" with respect to read counts. See
#' DESeq2 paper for dispersion model used.
#' @param logkdegsdtrend_slope Slope for log10(read count) vs. log(kdeg) replicate variability trend
#' @param logkdegsdtrend_intercept Intercept for log10(read count) vs. log(kdeg) replicate variability trend
#' @param logksynsdtrend_slope Slope for log10(read count) vs. log(ksyn) replicate variability trend
#' @param logksynsdtrend_intercept Intercept for log10(read count) vs. log(ksyn) replicate variability trend
#' @param logkdeg_sd Standard deviation of normal distribution from which reference log(kdeg)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logksyn_mean Mean of normal distribution from which reference log(ksyn)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logksyn_sd Standard deviation of normal distribution from which reference log(ksyn)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logkdeg_diff_avg Mean of normal distribution from which non-reference log(kdeg)
#' linear model parameters are drawn from for each feature if `param_details` is not provided.
#' @param logkdeg_diff_sd Standard deviation of normal distribution from which reference log(kdeg)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param logksyn_diff_avg Mean of normal distribution from which reference log(ksyn)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param logksyn_diff_sd Standard deviation of normal distribution from which reference log(ksyn)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param pdiff_kd Proportion of features for which non-reference log(kdeg) linear model parameters
#' differ from the reference.
#' @param pdiff_ks Proportion of features for which non-reference log(ksyn) linear model parameters
#' differ from the reference.
#' @param pdiff_both Proportion of features for which BOTH non-reference log(kdeg) and log(ksyn) linear model parameters
#' differ from the reference.
#' ksyns are simulated
#' @param pdo Dropout rate; think of this as the probability that a s4U containing
#' molecule is lost during library preparation and sequencing. If `pdo` is 0 (default)
#' then there is not dropout.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
EZSimulate <- function(nfeatures, ntreatments = 2, nreps = 3,
                       metadf = NULL,
                       mean_formula = NULL,
                       param_details = NULL,
                       seqdepth = 10000000, label_time = 2,
                       pnew = 0.05, pold = 0.001,
                       readlength = 200, Ucont = 0.25,
                       feature_prefix = "Gene",
                       dispslope = 5, dispint = 0.01,
                       logkdegsdtrend_slope = -0.3,
                       logkdegsdtrend_intercept = -2.25,
                       logksynsdtrend_slope = -0.3,
                       logksynsdtrend_intercept = -2.25,
                       logkdeg_mean = -1.9, logkdeg_sd = 0.7,
                       logksyn_mean = 2.3, logksyn_sd = 0.7,
                       logkdeg_diff_avg = 0, logksyn_diff_avg = 0,
                       logkdeg_diff_sd = 0.5, logksyn_diff_sd = 0.5,
                       pdiff_kd = 0.1, pdiff_ks = 0, pdiff_both = 0,
                       pdo = 0){


  ### NOTE: This is obviously currently a very trivial wrapper that could
  ### just be the default behavior of SimulateMultiCondition(). That being said,
  ### I suspect that there will be other aspects of simulation parameter setting
  ### that I would like to automate, so for now I will keep this as is

  ### Set parameters
  if(is.null(metadf)){

    mean_formula <- stats::as.formula('~treatment')
    metadf <- dplyr::tibble(sample = paste0('sample', 1:(nreps*ntreatments)),
                            treatment = rep(paste0('treatment', 1:ntreatments),
                                            each = nreps))


    simdata <- SimulateMultiCondition(nfeatures = nfeatures, metadf = metadf,
                                      mean_formula = mean_formula,
                                      param_details = param_details,
                                      seqdepth = seqdepth, label_time = label_time,
                                      pnew = pnew, pold = pold, readlength = readlength,
                                      Ucont = Ucont, feature_prefix = feature_prefix,
                                      dispslope = dispslope, dispint = dispint,
                                      logkdegsdtrend_slope = logkdegsdtrend_slope,
                                      logkdegsdtrend_intercept = logkdegsdtrend_intercept,
                                      logksynsdtrend_slope = logksynsdtrend_slope,
                                      logksynsdtrend_intercept = logksynsdtrend_intercept,
                                      logkdeg_mean = logkdeg_mean, logkdeg_sd = logkdeg_sd,
                                      logksyn_mean = logksyn_mean, logksyn_sd = logksyn_sd,
                                      logkdeg_diff_avg = logkdeg_diff_avg,
                                      logksyn_diff_avg = logksyn_diff_avg,
                                      logkdeg_diff_sd = logkdeg_diff_sd,
                                      logksyn_diff_sd = logksyn_diff_sd,
                                      pdiff_kd = pdiff_kd, pdiff_ks = pdiff_ks,
                                      pdiff_both = pdiff_both, pdo = pdo)


  }else{

    simdata <- SimulateMultiCondition(nfeatures = nfeatures, metadf = metadf,
                                      mean_formula = mean_formula,
                                      param_details = param_details,
                                      seqdepth = seqdepth, label_time = label_time,
                                      pnew = pnew, pold = pold, readlength = readlength,
                                      Ucont = Ucont, feature_prefix = feature_prefix,
                                      dispslope = dispslope, dispint = dispint,
                                      logkdegsdtrend_slope = logkdegsdtrend_slope,
                                      logkdegsdtrend_intercept = logkdegsdtrend_intercept,
                                      logksynsdtrend_slope = logksynsdtrend_slope,
                                      logksynsdtrend_intercept = logksynsdtrend_intercept,
                                      logkdeg_mean = logkdeg_mean, logkdeg_sd = logkdeg_sd,
                                      logksyn_mean = logksyn_mean, logksyn_sd = logksyn_sd,
                                      logkdeg_diff_avg = logkdeg_diff_avg,
                                      logksyn_diff_avg = logksyn_diff_avg,
                                      logkdeg_diff_sd = logkdeg_diff_sd,
                                      logksyn_diff_sd = logksyn_diff_sd,
                                      pdiff_kd = pdiff_kd, pdiff_ks = pdiff_ks,
                                      pdiff_both = pdiff_both, pdo = pdo)


  }



}





#' Simulate NR-seq data for multiple replicates of multiple biological conditions
#'
#' `SimulateMultiCondition` is a highly flexibly simulator that combines linear modeling
#' of log(kdeg)'s and log(ksyn)'s with `SimulateOneRep` to simulate an NR-seq dataset. The linear model
#' allows you to simulate multiple distinct treatments, batch effects, interaction effects,
#' etc. The current downside for its flexibility is its relative complexity to implement.
#' Easier to use simulators are on the way to EZbakR.
#'
#' @param nfeatures Number of "features" (e.g., genes) to simulate data for
#' @param metadf A data frame with the following columns:
#' \itemize{
#'  \item sample: Names given to samples to simulate.
#'  \item <details>: Any number of columns with any names (not taken by other metadf columns)
#'  storing factors by which the samples can be stratified. These can be referenced
#'  in `mean_formula`, described below.
#' }
#' These parameters (described more below) can also be included in metadf to specify sample-specific simulation
#' parameter:
#' \itemize{
#'  \item seqdepth
#'  \item label_time
#'  \item pnew
#'  \item pold
#'  \item readlength
#'  \item Ucont
#' }
#' @param mean_formula A formula object that specifies the linear model used to
#' relate the factors in the <details> columns of `metadf` to average log(kdegs) and
#' log(ksyns) in each sample.
#' @param param_details A data frame with one row for each column of the design matrix
#' obtained from `model.matrix(mean_formula, metadf)` that describes how to simulate
#' the linear model parameters. The columns of this data frame are:
#' \itemize{
#'  \item param: Name of linear model parameter as it appears in the column names of the
#'  design matrix from `model.matrix(mean_formula, metadf)`.
#'  \item reference: Boolean; TRUE if you want to treat that parameter as a "reference". This
#'  means that all other parameter values that aren't global parameters are set equal to this
#'  unless otherwise determined (see `pdiff_*` parameters for how it is determined if a parameter
#'  will differ from the reference).
#'  \item global: Boolean; TRUE if you want to treat that parameter as a global parameter. This means
#'  that a single value is used for all features.
#'  \item logkdeg_mean: If parameter is the reference, then its value for the log(kdeg) linear model
#'  will be drawn from a normal distribution with this mean. If it is a global parameter, then this
#'  value will be used. If it is neither of these, then its value in the log(kdeg) linear model will
#'  either be the reference (if there is no difference between this condition's value and the reference)
#'  or the reference's value + a normally distributed random variable centered on this value.
#'  \item logkdeg_sd: sd used for draws from normal distribution as described for `logkdeg_mean`.
#'  \item logksyn_mean: Same as `logkdeg_mean` but for log(ksyn) linear model.
#'  \item logksyn_sd: Same as `logkdeg_sd` but for log(kdeg) linear model.
#'  \item pdiff_ks: Proportion of features whose value of this parameter in the log(ksyn) linear model
#'  will differ from the reference's. Should be a number between 0 and 1, inclusive. For example, if
#'  `pdiff_ks` is 0.1, then for 10% of features, this parameter will equal the reference parameter +
#'  a normally distributed random variable with mean `logksyn_mean` and sd `logksyn_sd`. For the other
#'  90% of features, this parameter will equal the reference.
#'  \item pdiff_kd: Same as `pdiff_ks` but for log(kdeg) linear model.
#'  \item pdiff_both: Proportion of features whose value for this parameter in BOTH the
#'  log(kdeg) and log(ksyn) linear models will differ from the reference. Value must be
#'  between 0 and min(c(pdiff_kd, pdiff_ks)) in that row.
#' }
#' If param_details is not specified by the user, the first column of the design matrix
#' is assumed to represent the reference parameter, all parameters are assumed to be
#' non-global, logkdeg_mean and logksyn_mean are set to the equivalently named parameter values
#' described below for the reference and `logkdeg_diff_avg` and `logksyn_diff_avg` for all other parameters,
#' logkdeg_sd and logksyn_sd are set to the equivalently named parameter values
#' described below for the reference and `logkdeg_diff_sd` and `logksyn_diff_sd` for all other parameters,
#' and pdiff_kd, pdiff_ks, and pdiff_both are all set to the equivalently named parameter values.
#' @param seqdepth Only relevant if `read_vect` is not provided; in that case, this is
#' the total number of reads to simulate.
#' @param label_time Length of s^{4}U feed to simulate.
#' @param pnew Probability that a T is mutated to a C if a read is new.
#' @param pold Probability that a T is mutated to a C if a read is old.
#' @param readlength Length of simulated reads. In this simple simulation, all reads
#' are simulated as being exactly this length.
#' @param Ucont Probability that a nucleotide in a simulated read is a U.
#' @param logkdeg_mean Mean of normal distribution from which reference log(kdeg)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param feature_prefix Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows up in the
#' `feature` column of the output simulated data table.
#' @param dispslope Negative binomial dispersion parameter "slope" with respect to read counts. See
#' DESeq2 paper for dispersion model used.
#' @param dispint Negative binomial dispersion parameter "intercept" with respect to read counts. See
#' DESeq2 paper for dispersion model used.
#' @param logkdegsdtrend_slope Slope for log10(read count) vs. log(kdeg) replicate variability trend
#' @param logkdegsdtrend_intercept Intercept for log10(read count) vs. log(kdeg) replicate variability trend
#' @param logksynsdtrend_slope Slope for log10(read count) vs. log(ksyn) replicate variability trend
#' @param logksynsdtrend_intercept Intercept for log10(read count) vs. log(ksyn) replicate variability trend
#' @param logkdeg_sd Standard deviation of normal distribution from which reference log(kdeg)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logksyn_mean Mean of normal distribution from which reference log(ksyn)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logksyn_sd Standard deviation of normal distribution from which reference log(ksyn)
#' linear model parameter is drawn from for each feature if `param_details` is not provided.
#' @param logkdeg_diff_avg Mean of normal distribution from which non-reference log(kdeg)
#' linear model parameters are drawn from for each feature if `param_details` is not provided.
#' @param logkdeg_diff_sd Standard deviation of normal distribution from which reference log(kdeg)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param logksyn_diff_avg Mean of normal distribution from which reference log(ksyn)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param logksyn_diff_sd Standard deviation of normal distribution from which reference log(ksyn)
#' linear model parameter are drawn from for each feature if `param_details` is not provided.
#' @param pdiff_kd Proportion of features for which non-reference log(kdeg) linear model parameters
#' differ from the reference.
#' @param pdiff_ks Proportion of features for which non-reference log(ksyn) linear model parameters
#' differ from the reference.
#' @param pdiff_both Proportion of features for which BOTH non-reference log(kdeg) and log(ksyn) linear model parameters
#' differ from the reference.
#' ksyns are simulated
#' @param pdo Dropout rate; think of this as the probability that a s4U containing
#' molecule is lost during library preparation and sequencing. If `pdo` is 0 (default)
#' then there is not dropout.
#' @import data.table
#' @importFrom magrittr %>%
#' @export
SimulateMultiCondition <- function(nfeatures, metadf, mean_formula,
                                   param_details = NULL,
                                   seqdepth = 10000000, label_time = 2,
                                   pnew = 0.05, pold = 0.001,
                                   readlength = 200, Ucont = 0.25,
                                   feature_prefix = "Gene",
                                   dispslope = 5, dispint = 0.01,
                                   logkdegsdtrend_slope = -0.3,
                                   logkdegsdtrend_intercept = -2.25,
                                   logksynsdtrend_slope = -0.3,
                                   logksynsdtrend_intercept = -2.25,
                                   logkdeg_mean = -1.9, logkdeg_sd = 0.7,
                                   logksyn_mean = 2.3, logksyn_sd = 0.7,
                                   logkdeg_diff_avg = 0, logksyn_diff_avg = 0,
                                   logkdeg_diff_sd = 0.5, logksyn_diff_sd = 0.5,
                                   pdiff_kd = 0.1, pdiff_ks = 0, pdiff_both = 0,
                                   pdo = 0){

  ### Hack to deal with devtools::check() NOTEs
  reference <- param <- NULL


  `.` <- list



  ### Create param_details if not provided

  mean_design <- stats::model.matrix(mean_formula, metadf)

  mean_design_cols <- colnames(mean_design)


  if(is.null(param_details)){

    lmdc <- length(mean_design_cols)

    param_details <- dplyr::tibble(
      param = mean_design_cols,
      reference = rep(c(TRUE, FALSE), times = c(1, lmdc - 1)),
      global = FALSE,
      logkdeg_mean = rep(c(logkdeg_mean, logkdeg_diff_avg), times = c(1, lmdc - 1)),
      logkdeg_sd = rep(c(logkdeg_sd, logkdeg_diff_sd), times = c(1, lmdc - 1)),
      logksyn_mean = rep(c(logksyn_mean, logksyn_diff_avg), times = c(1, lmdc - 1)),
      logksyn_sd = rep(c(logksyn_sd, logksyn_diff_sd), times = c(1, lmdc - 1)),
      pdiff_kd = pdiff_kd,
      pdiff_ks = pdiff_ks,
      pdiff_both = pdiff_both
    )

  }



  ### Check validity of input

  args <- c(as.list(environment()))

  check_SimulateMultiCondition_input(args)


  ### Fill metadf with parameters that are only specified as single value

  mcols <- colnames(metadf)

  if(!('seqdepth' %in% mcols)){

    metadf$seqdepth <- seqdepth

  }

  if(!('label_time' %in% mcols)){

    metadf$label_time <- label_time

  }

  if(!('pnew' %in% mcols)){

    metadf$pnew <- pnew

  }


  if(!('pold' %in% mcols)){

    metadf$pold <- pold

  }

  if(!('pdo' %in% mcols)){

    metadf$pdo <- pdo

  }






  ### Need to simulate linear model parameter values for all parameters specified


  # Reference log(kdegs)
  pdref <- param_details %>%
    dplyr::filter(reference)

  logkdeg_ref <- stats::rnorm(nfeatures,
                       mean = pdref$logkdeg_mean,
                       sd = pdref$logkdeg_sd)


  # Reference log(ksyns)
  logksyn_ref <- stats::rnorm(nfeatures,
                       mean = pdref$logksyn_mean,
                       sd = pdref$logksyn_sd)


  # Number of samples to simulate; to be used multiple times later
  nsamp <- nrow(metadf)

  # Preallocate lists to store log(kdeg) and log(ksyn)s linear model parameters
  # for each sample and feature
  logkdeg_params <- vector(mode = "list",
                           length = length(mean_design_cols))
  logksyn_params <- logkdeg_params

  # Determine linear model parameter values
  for(p in seq_along(mean_design_cols)){

    pd <- param_details %>%
      dplyr::filter(param == mean_design_cols[p])

    if(pd$reference){

      ### EASY: Just use reference value for reference parameters
      logkdeg_params[[p]] <- logkdeg_ref
      logksyn_params[[p]] <- logksyn_ref

    }else if(pd$global){

      ### EASY: Just use global value for global parameters
      logkdeg_params[[p]] <- rep(pd$logkdeg_mean, times = nfeatures)
      logksyn_params[[p]] <- rep(pd$logksyn_mean, times = nfeatures)

    }else{

      ### HARD: Need to simulate non-reference, non-global values with respect
      ### to the references. Some fraction of these will differ from reference,
      ### some fraction will be exactly the same as the reference.
      ###
      ### User also specifies what fraction of the time they want both kdeg and
      ### ksyn to differ with respect to reference.

      ndiff_kd <- pd$pdiff_kd*nfeatures
      ndiff_ks <- pd$pdiff_ks*nfeatures
      ndiff_both <- pd$pdiff_both*nfeatures

      diff_kd_end <- ceiling(ndiff_kd)
      diff_ks_start <- ceiling(diff_kd_end - ndiff_both)
      diff_ks_end <- diff_ks_start + ndiff_ks

      is_kdeg_param_nonzero <- rep(c(1, 0),
                                   times = c(diff_kd_end,
                                             nfeatures - diff_kd_end))

      is_ksyn_param_nonzero <- rep(c(0, 1, 0),
                                   times = c(diff_ks_start,
                                             ndiff_ks,
                                             nfeatures - diff_ks_start - ndiff_ks))

      logkdeg_params[[p]] <- is_kdeg_param_nonzero*stats::rnorm(nfeatures,
                                                         pd$logkdeg_mean,
                                                         pd$logkdeg_sd) +
        logkdeg_ref

      logksyn_params[[p]] <- is_ksyn_param_nonzero*stats::rnorm(nfeatures,
                                                         pd$logksyn_mean,
                                                         pd$logksyn_sd) +
        logksyn_ref


    }

  }

  # Function to extract the ith element of a vector
  # For grabbing a feature's set of parameter from logkdeg/ksyn_params
  extract_ith <- function(list, i) {
    sapply(list, function(x) x[i])
  }

  # Function to compute mean log(kdeg), log(ksyn), and abundance
  # in each sample.
  compute_kinetics <- function(X, logkdeg_params,
                               logksyn_params,
                               pdo,
                               tl,
                               n = nfeatures) {

    abundance <- vector("list", n)
    logkdeg <- abundance
    logksyn <- abundance

    for (i in 1:n) {
      logksyn_i <- X %*% extract_ith(logksyn_params, i)
      logkdeg_i <- X %*% extract_ith(logkdeg_params, i)

      logkdeg[[i]] <- logkdeg_i
      logksyn[[i]] <- logksyn_i

      # Compute dropout adjusted abundance; just ksyn/kdeg if no dropout
      fn <- 1 - exp(-exp(logkdeg_i)*tl)
      abundance[[i]] <- (exp(logksyn_i)/exp(logkdeg_i))*((1 - fn) + (1 - pdo)*fn)

    }

    return(list(abundance = abundance,
                logkdeg = logkdeg,
                logksyn = logksyn))
  }

  # Function to normalize abundances and calculate expected read count
  # for each feature in each sample
  compute_readcounts <- function(abundances,
                                 seqdepths,
                                 n = nsamp){


    sums <- rep(0, times = n)
    means <- sums
    sds <- sums

    for(i in 1:n){

      abundances_i <- extract_ith(abundances, i)



      sums[i] <- sum(abundances_i)
      means[i] <- mean(log10(abundances_i))
      sds[i] <- stats::sd(log10(abundances_i))

    }

    readcounts <- lapply(abundances,
                         function(vector) (vector/sums)*seqdepths)
    zscores <- lapply(abundances,
                      function(vector) (log10(vector) - means)/sds)
    return(list(readcounts = readcounts,
                read_zscores = zscores))


  }


  ### Determine sample-specific averages to simulate from

  kinetics <- compute_kinetics(mean_design, logkdeg_params,
                               logksyn_params,
                               pdo = metadf$pdo,
                               tl = metadf$label_time)


  reads <- compute_readcounts(kinetics$abundance,
                              seqdepths = metadf$seqdepth)


  kinetics_and_reads <- c(kinetics, reads)


  ### Simulate per-replicate values using calculated means
  # log(kdeg) ~ Normal()
  # log(ksyn) ~ Normal()
  # read count ~ Negative Binomial()

  logkdegs <- vector(mode = "list",
                     length = nfeatures)
  logksyns <- logkdegs
  reads <- logkdegs

  for(i in 1:nfeatures){

    feature_logkdegs <- kinetics_and_reads$logkdeg[[i]]
    feature_logksyns <- kinetics_and_reads$logksyn[[i]]
    feature_readavgs <- kinetics_and_reads$readcounts[[i]]

    logkdeg_sds <- exp(kinetics_and_reads$read_zscores[[i]]*logkdegsdtrend_slope +
                         logkdegsdtrend_intercept)

    logksyn_sds <- exp(kinetics_and_reads$read_zscores[[i]]*logksynsdtrend_slope +
                         logksynsdtrend_intercept)



    ### Each element is a vector of sample-specific kinetic parameters and read counts
    # NOTE: bit weird to define sample-specific kinetic parameters rather than
    # sample specific fraction news. From a practical perspective though, this
    # is equivalent to that strategy.
    logkdegs[[i]] <- stats::rnorm(n = nsamp,
                           mean = feature_logkdegs,
                           sd = logkdeg_sds)

    logksyns[[i]] <- stats::rnorm(n = nsamp,
                           mean = feature_logksyns,
                           sd = logksyn_sds)

    reads[[i]] <- stats::rnbinom(n = nsamp,
                          mu = feature_readavgs,
                          size = 1/(dispslope/feature_readavgs + dispint))


  }


  ### Simulate data for each replicate
  simdata <- vector(mode = "list", length = nrow(metadf))
  for(s in 1:nrow(metadf)){


    pdo <- metadf$pdo[s]

    kdeg_vect <- exp(extract_ith(logkdegs, s))
    fn_vect <- 1 - exp(-kdeg_vect*metadf$label_time[s])

    # Dropout adjusted
    fn_vect <- (fn_vect*(1 - pdo))/(fn_vect*(1 - pdo) + (1 - fn_vect))
    kdeg_vect <- -log(1 - fn_vect)/metadf$label_time[s]



    # A bit unclear just from the interface (so suboptimal function design)
    # but the key here is that if there is dropout, fn_vect will represent
    # the dropout biased fraction new, and kdeg_vect will represent the
    # kdeg that would be estimated from the true, unbiased fraction new.
    # This ensures that the ground_truth table from this function contains
    # both the biased and unbiased fraction new estimates in it.
    simdata[[s]] <- SimulateOneRep(nfeatures = nfeatures,
                                   read_vect = extract_ith(reads, s),
                                   label_time = as.numeric(metadf[s,"label_time"]),
                                   sample_name = as.character(metadf[s, "sample"]),
                                   fn_vect = fn_vect,
                                   kdeg_vect = exp(extract_ith(logkdegs, s)),
                                   ksyn_vect = exp(extract_ith(logksyns, s)),
                                   pnew = as.numeric(metadf[s, "pnew"]),
                                   pold = as.numeric(metadf[s, "pold"]))

  }

  # Combine replicate simulation data objects into one
  names_to_bind <- names(simdata[[1]])

  final_simdata <- lapply(names_to_bind, function(name) {
    dplyr::bind_rows(lapply(simdata, function(inner_list) inner_list[[name]]))
  })

  names(final_simdata) <- names_to_bind


  ### Save linear model parameters in convenient format

  names(logkdeg_params) <- paste0("true_logkdeg_", mean_design_cols)
  names(logksyn_params) <- paste0("true_logksyn_", mean_design_cols)

  feature_prefix <- "Gene"

  kinetic_parameters <- dplyr::bind_cols(list(logkdeg_params, logksyn_params))
  kinetic_parameters[['feature']] <- paste0(feature_prefix, 1:nfeatures)



  ### Gather output

  # TO-DO: Am not tracking the replicate fraction news in a way
  # that accounts for potential dropout. PerRepTruth will only
  # contain the dropout biased fraction new estimate
  if(all(metadf$pdo == 0)){

    output <- list(cB = final_simdata$cB,
                   PerRepTruth = final_simdata$ground_truth,
                   AvgTruth = kinetic_parameters,
                   metadf = metadf,
                   param_details = param_details)

  }else{

    nodropout_truth <- vector(mode = "list", length = nrow(metadf))
    for(s in 1:nrow(metadf)){

      pdo <- metadf$pdo[s]

      kdeg_vect <- exp(extract_ith(logkdegs, s))
      fn_vect <- 1 - exp(-kdeg_vect*metadf$label_time[s])


      nodropout_truth[[s]] <- data.table::data.table(sample = metadf$sample[s],
                                      feature = paste0(feature_prefix, 1:nfeatures),
                                      unbiased_fraction_highTC = fn_vect,
                                      unbiased_kdeg = kdeg_vect)

    }

    output <- list(cB = final_simdata$cB,
                   PerRepTruth = final_simdata$ground_truth,
                   UnbiasedFractions = dplyr::bind_rows(nodropout_truth),
                   AvgTruth = kinetic_parameters,
                   metadf = metadf,
                   param_details = param_details)

  }


  return(output)


}



#' Simulation of transcript isoform kinetic parameters.
#'
#' `SimulateIsoforms()` performs a simple simulation of isoform-specific kinetic
#' parameters to showcase and test `EstimateIsoformFractions()`. It assumes that
#' there are a set of reads (fraction of total set by `funique` parameter) which
#' map uniquely to a given isoform, while the rest are ambiguous to all isoforms
#' from that gene. Mutational content of these reads are simulated as in
#' `SimulateOneRep()`.
#'
#' @param nfeatures Number of "features" to simulate data for. Each feature will
#' have a simulated number of transcript isoforms
#' @param nt (Optional), can provide a vector of the number of isoforms you would
#' like to simulate for each of the `nfeatures` features. Vector can either be length
#' 1, in which case that many isoforms will be simulated for all features, or length
#' equal to `nfeatures`.
#' @param seqdepth Total number of sequencing reads to simulate
#' @param label_time Length of s^{4}U feed to simulate.
#' @param feature_prefix Name given to the i-th feature is `paste0(feature_prefix, i)`. Shows up in the
#' `feature` column of the output simulated data table.
#' @param sample_name Character vector to assign to `sample` column of output simulated
#' data table (the cB table).
#' @param pnew Probability that a T is mutated to a C if a read is new.
#' @param pold Probability that a T is mutated to a C if a read is old.
#' @param readlength Length of simulated reads. In this simple simulation, all reads
#' are simulated as being exactly this length.
#' @param Ucont Probability that a nucleotide in a simulated read is a U.
#' @param avg_numiso Average number of isoforms for each feature. Feature-specific
#' isoform counts are drawn from a Poisson distribution with this average. NOTE:
#' to insure that all features have multiple isoforms, the simulated number of
#' isoforms drawn from a Poisson distribution is incremented by 2. Thus, the
#' actual average number of isoforms from each feature is `avg_numiso` + 2.
#' @param psynthdiff Percentage of genes for which all isoform abundance differences
#' are synthesis driven. If not synthesis driven, then isoform abundance differences
#' will be driven by differences in isoform kdegs.
#' @param logkdeg_mean meanlog of a log-normal distribution from which
#' kdegs are simulated
#' @param logkdeg_sd sdlog of a log-normal distribution from which
#' kdegs are simulated
#' @param logksyn_mean meanlog of a log-normal distribution from which
#' ksyns are simulated
#' @param logksyn_sd sdlog of a log-normal distribution from which
#' ksyns are simulated
#' @export
SimulateIsoforms <- function(nfeatures,
                             nt = NULL,
                             seqdepth = nfeatures*2500,
                             label_time = 4,
                             sample_name = 'sampleA',
                             feature_prefix = 'Gene',
                             pnew = 0.1,
                             pold = 0.002,
                             funique = 0.2,
                             readlength = 200,
                             Ucont = 0.25,
                             avg_numiso = 2,
                             psynthdiff = 0.5,
                             logkdeg_mean = -1.9, logkdeg_sd = 0.7,
                             logksyn_mean = 2.3, logksyn_sd = 0.7
                             ){


  if(is.null(nt)){

    # Number of isoforms per genes (make them all multi-isoformgenes)
    nt <- stats::rpois(nfeatures, avg_numiso) + 2

  }else if(length(nt) == 1){

    nt <- rep(nt, times = nfeatures)

  }

  # Is synthesis different
  syn_driven <- as.logical(stats::rbinom(nfeatures, size = 1, prob = psynthdiff))

  # Proportions for each transcript
  pt <- rep(0, times = sum(nt))

  tracker <- 1
  for(i in 1:nfeatures){

    # Simulate one dominant isoform and then a number of less prevalent isoforms
    # Idea is to generate a number on an unnormalized scale (9 = dominant,
    # other isoforms randomly assigned a real number between 1 and 4, uniform
    # distribution over that range), and then normalize it for each gene (i.e,
    # divide unnormalized number by sum of unnormalized numbers for that gene)
    pt[tracker:(tracker+(nt[i] - 1))] <- c(9, stats::runif(n = (nt[i] - 1), min = 1, max = 4))
    pt[tracker:(tracker+(nt[i] - 1))] <- pt[tracker:(tracker+(nt[i] - 1))]/sum(pt[tracker:(tracker+(nt[i] - 1))])

    tracker <- tracker + nt[i]
  }


  # ksyn for each gene
  ks_g <- stats::rlnorm(nfeatures, meanlog = logksyn_mean,
                 sdlog = logksyn_sd)

  # kdeg for each gene
  kd_g <- stats::rlnorm(nfeatures, meanlog = logkdeg_mean,
                 sdlog = logkdeg_sd)

  # Simulate ks and kdeg for each transcript -------------------------------------

  gene_vect <- rep(1:nfeatures, times = nt)

  ks_t <- rep(0, times = sum(nt))
  kd_t <- rep(0, times = sum(nt))

  for(i in 1:sum(nt)){

    dom_p <- max(pt[gene_vect == gene_vect[i]])

    # Dominant isoform is given the simulated gene-wide kinetic parameters
    # Each non-dominant isoform is either termed a "synthesis-driven" alternative
    # isoform or a "degradation-driven". A synthesis-driven alt. isoform is one
    # whose lower than dominant abundance is driven by a lower rate of synthesis
    # for that gene. A degradation-driven alt. isoform is one whose lower than
    # dominant abundance is driven by lower stability. Synthesis and degradation
    # rate constants are chosen so that steady-state abundances are consistent
    # with the isoform percentages simulated earlier
    if(pt[i] == dom_p){
      ks_t[i] <- ks_g[gene_vect[i]]
      kd_t[i] <- kd_g[gene_vect[i]]
    }else{
      if(syn_driven[gene_vect[i]]){
        kd_t[i] <- kd_g[gene_vect[i]]

        ks_t[i] <- ((pt[i])/dom_p)*ks_g[gene_vect[i]]
      }else{

        ks_t[i] <- ks_g[gene_vect[i]]

        kd_t[i] <- (dom_p/pt[i])*kd_g[gene_vect[i]]

      }

    }

  }


  # Simulate reads from each transcript ------------------------------------------

  ### Non-unique reads
  # I randomly choose some fraction of reads to be uniquely mapping to one isoform
  # and some fraction to map to multiple isoforms.
  # For the non-unique isoforms, I assume that the probability a read mapped to
  # a given transcript is simply that transcript's relative abundance (i.e.,
  # the proportion of RNA generated from that gene which become that isoform.
  # this is `pt` simulated earlier).
  # NOTE:
  # I got confused once about how the mutational data is simulated in the case
  # of non-unique reads. These reads are unambiguously from a given transcript
  # isoform, so the fraction new that drives the simulation of the mutational
  # data is just that isoform's simulated fraction new. At one point I wondered
  # why I wasn't probablistically determining the transcript of origin for
  # a non-unique read. This is unnecessary as the number of reads originating
  # from a given isoform takes into account the relative abundance of each isoform.

  if(funique < 1){
    reads <- ceiling(((ks_t/kd_t)/sum(ks_t/kd_t))*seqdepth*(1-funique))
    fns <- rep(1 - exp(-kd_t*label_time), times = reads)


    # Number of Us
    nUs <- stats::rbinom(n = sum(reads), size = readlength, prob = Ucont)

    # Newness
    newness <- stats::rbinom(n = sum(reads), size = 1, p = fns)

    # Mutations
    TCs <- stats::rbinom(n = sum(reads), size = nUs,
                         prob = newness*pnew + (1 - newness)*pold )

    # Transcript IDs
    transcript_vect <- rep('', times = nfeatures)
    for(i in 1:nfeatures){
      ### TO-DO: Make this a prefix a user provided parameter
      transcript_set <- paste0(feature_prefix, i, "_", "Transcript", 1:nt[i])
      transcript_vect[i] <- paste(transcript_set, collapse = "+")

    }

    # Make cB
    cB <- dplyr::tibble(TC = TCs, nT = nUs, read_ID = 1:length(TCs),
                 GF = paste0(feature_prefix, rep(gene_vect, times = reads)),
                 transcripts = rep(rep(transcript_vect, times = nt), times = reads))

  }else{
    cB <- dplyr::tibble()
  }



  ### Unique reads


  reads_u <- ceiling(((ks_t/kd_t)/sum(ks_t/kd_t))*seqdepth*(funique))
  fns <- rep(1 - exp(-kd_t*label_time), times = reads_u)


  # Number of Us
  nUs <- stats::rbinom(n = sum(reads_u), size = readlength, prob = Ucont)

  # Newness
  newness <- stats::rbinom(n = sum(reads_u), size = 1, p = fns)

  # Mutations
  TCs <- stats::rbinom(n = sum(reads_u), size = nUs, prob = newness*pnew + pold - newness*pold )

  # Transcript IDs
  gene_vect <- rep(1:nfeatures, times = nt)
  transcript_vect <- c()
  for(i in 1:nfeatures){
    ### TO-DO: Make this a prefix a user provided parameter
    transcript_set <- paste0(feature_prefix, i, "_", "Transcript", 1:nt[i])
    transcript_vect <- c(transcript_vect, transcript_set)

  }

  # Make cB
  cB_u <- tibble(TC = TCs, nT = nUs, read_ID = (nrow(cB) + 1):(nrow(cB) + length(TCs)),
                 GF = paste0(feature_prefix, rep(gene_vect, times = reads_u)),
                 transcripts = rep(transcript_vect, times = reads_u))


  if(funique >= 1){
    # Merge
    cB <- cB_u
  }else{
    # Merge
    cB <- dplyr::bind_rows(cB, cB_u)
  }


  ### Assemble ground truth and data

  cB <- cB %>%
    dplyr::rename(feature = GF) %>%
    dplyr::group_by(feature, transcripts, TC, nT) %>%
    dplyr::count() %>%
    dplyr::mutate(sample = sample_name) %>%
    dplyr::select(sample, feature, transcripts, TC, nT, n)

  truth <- dplyr::tibble(
    feature = paste0(feature_prefix, gene_vect),
    transcript_id = transcript_vect,
    true_kdeg = kd_t,
    true_ksyn = ks_t,
    true_fn = 1 - exp(-kd_t*label_time),
    true_count = seqdepth*ceiling((ks_t/kd_t)/sum(ks_t/kd_t)),
    true_TPM = (ks_t/kd_t)/(sum(ks_t/kd_t)/1000000)
  )

  return(list(cB = cB,
              ground_truth = truth))


}



#########################
###### PARAMETER CHECKS #
#########################


check_SimulateOneRep_input <- function(args){

  ### nfeatures

  NF <- args$nfeatures

  if(!is.numeric(NF)){

    stop("nfeatures must be numeric!")

  }

  if(NF < 1){
    stop("nfeatures must be >= 1!")
  }

  if(round(NF) != NF){

    stop("nfeatures must be an integer!")

  }


  ### read_vect
  rv <- args$read_vect

  if(!is.null(rv)){

    if(length(rv) != 1 & length(rv) != NF){

      stop("read_vect must be either length 1 or length nfeatures!")

    }

    if(!all(is.numeric(rv))){

      stop("All elements of read_vect must be numeric!")

    }

    if(!all(rv >= 0)){

      stop("All elements of read_vect must be >= 0!")
    }

    if(!all(round(rv) == rv)){

      stop("All elements of read_vect must be integers!")

    }

  }


  ### label_time
  tl <- args$label_time

  if(!is.numeric(tl)){

    stop("label_time must be numeric")

  }

  if(tl < 0){

    stop("label_time must be >= 0!")

  }


  ### sample_name
  sname <- args$sample_name

  if(!is.character(sname)){

    stop("sample_name should be a string!")

  }

  if(length(sname) > 1){

    stop("sample_name should be a single string!")

  }

  ### feature_prefix
  fp <- args$feature_prefix

  if(!is.character(fp)){

    stop("feature_prefix should be a string!")

  }

  if(length(fp) > 1){

    stop("feature_prefix should be a single string!")

  }

  ### logkdeg_mean
  lkd <- args$logkdeg_mean

  if(!is.numeric(lkd)){

    stop("logkdeg_mean must be numeric")

  }


  ### logkdeg_sd
  lkd_sd <- args$logkdeg_sd

  if(!is.numeric(lkd_sd)){

    stop("logkdeg_mean must be numeric!")

  }

  if(lkd_sd <= 0){

    stop("logkdeg_sd must be >= 0")

  }


  ### logksyn_mean
  lks <- args$logksyn_mean

  if(!is.numeric(lks)){

    stop("logksyn_mean must be numeric!")

  }


  ### logkdeg_sd
  lks_sd <- args$logksyn_sd

  if(!is.numeric(lks_sd)){

    stop("logksyn_mean must be numeric!")

  }

  if(lks_sd <= 0){

    stop("logksyn_sd must be >= 0")

  }


  ### seqdepth
  sdep <- args$seqdepth

  if(!is.numeric(sdep)){

    stop("seqdepth must be numeric")

  }

  if(sdep <= 0){

    stop("seqdepth must be > 0")

  }

  if(round(sdep) != sdep){

    stop("seqdepth must be an integer!")

  }


  ### pnew and old
  pnew <- args$pnew
  pold <- args$pold

  if(!is.numeric(pnew)){

    stop("pnew must be numeric!")

  }

  if(!is.numeric(pold)){

    stop("pold must be numeric!")

  }

  if(pnew <= 0){

    stop("pnew must be > 0")

  }

  if(pold < 0){

    stop("pnew must be >= 0")

  }

  if(pnew > 1){
    stop("pnew must be <= 1")
  }

  if(pold >= 1){
    stop("pold must be < 1")
  }

  if(pnew <= pold){
    stop("pnew must be strictly greater than pnew!")
  }



  ### Ucont
  Ucont <- args$Ucont

  if(!is.numeric(Ucont)){

    stop("Ucont must be numeric!")

  }

  if(Ucont <= 0){

    stop("Ucont must be > 0")

  }

  if(Ucont > 1){

    stop("Ucont must be <= 1")

  }


}



check_SimulateMultiCondition_input <- function(args){

  metadf <- args$metadf
  pd <- args$param_details

  ### nfeatures

  NF <- args$nfeatures

  if(!is.numeric(NF)){

    stop("nfeatures must be numeric!")

  }

  if(NF < 1){
    stop("nfeatures must be >= 1!")
  }

  if(round(NF) != NF){

    stop("nfeatures must be an integer!")

  }


  ### label_time
  tl <- args$label_time

  if(!is.numeric(tl)){

    stop("label_time must be numeric")

  }

  if(tl < 0){

    stop("label_time must be >= 0!")

  }


  ### sample_name
  sname <- args$metadf$sample

  if(!all(is.character(sname))){

    stop("All elements of metadf column `sample` must be strings!")

  }


  ### feature_prefix
  fp <- args$feature_prefix

  if(!is.character(fp)){

    stop("feature_prefix should be a string!")

  }

  if(length(fp) > 1){

    stop("feature_prefix should be a single string!")

  }


  ### logkdeg_sd
  lkd_sd <- args$logkdeg_sd

  if(!is.numeric(lkd_sd)){

    stop("logkdeg_mean must be numeric!")

  }

  if(lkd_sd <= 0){

    stop("logkdeg_sd must be >= 0")

  }


  ### seqdepth
  sdep <- args$seqdepth

  if(!is.numeric(sdep)){

    stop("seqdepth must be numeric")

  }

  if(sdep <= 0){

    stop("seqdepth must be > 0")

  }

  if(round(sdep) != sdep){

    stop("seqdepth must be an integer!")

  }


  ### pnew and old
  pnew <- args$pnew
  pold <- args$pold

  if(!is.numeric(pnew)){

    stop("pnew must be numeric!")

  }

  if(!is.numeric(pold)){

    stop("pold must be numeric!")

  }

  if(pnew <= 0){

    stop("pnew must be > 0")

  }

  if(pold < 0){

    stop("pnew must be >= 0")

  }

  if(pnew > 1){
    stop("pnew must be <= 1")
  }

  if(pold >= 1){
    stop("pold must be < 1")
  }

  if(pnew <= pold){
    stop("pnew must be strictly greater than pnew!")
  }



  ### Ucont
  Ucont <- args$Ucont

  if(!is.numeric(Ucont)){

    stop("Ucont must be numeric!")

  }

  if(Ucont <= 0){

    stop("Ucont must be > 0")

  }

  if(Ucont > 1){

    stop("Ucont must be <= 1")

  }


}
