#' Run quality control checks
#'
#' `EZQC()` assesses multiple aspects of your NR-seq data and generates a number
#' of plots visualizing dataset-wide trends.
#'
#' `EZQC()` checks the following aspects of your NR-seq data. If you have passed
#' an `EZbakRData` object, then `EZQC()` checks:
#' \itemize{
#'  \item Raw mutation rates: In all sequencing reads, how many T's in the reference
#'  were a C in the read? The hope is that raw mutation rates are higher than -label
#'  controls in all +label samples. Higher raw mutation rates, especially when using
#'  standard label times (e.g., 2 hours or more in mammalian systems), are typically
#'  a sign of good label incorporation and low labeled RNA/read dropout. If you
#'  don't have -label samples, know that background mutation rates are typically
#'  less than 0.2%, so +label raw mutation rates several times higher than this
#'  would be preferable.
#'  \item Mutation rates in labeled and unlabeled reads: The raw mutation rate
#'  counts all mutations in all reads. In a standard NR-seq experiment performed
#'  with a single metabolic label, there are typically two populations of reads:
#'  1) Those from labeled RNA, having higher mutation rates due to chemical conversion/recoding
#'  of the metabolic label and 2) those from unlabeled RNA, having lower, background
#'  levels of mutations. `EZbakR` fits a two component mixture model to estimate the mutation
#'  rates in these two populations separately. A successful NR-seq experiment should
#'  have a labeled read mutation rate of > 1% and a low background mutation rate of
#'  < 0.3%.
#'  \item Read count replicate correlation: Simply the log10 read count correlation
#'  for replicates, as inferred from your `metadf`.
#' }
#'
#' If you have passed an `EZbakRFractions` object, i.e., the output of `EstimateFractions()`,
#' then in addition to the checks in the `EZbakRData` input case, `EZQC()` also checks:
#' \itemize{
#'  \item Fraction labeled distribution: This is the distribution of feature-wise
#'  fraction labeled's (or fraction high mutation content's) estimated by `EstimateFractions()`.
#'  The "ideal" is a distribution with mean around 0.5, as this maximizes the amount of RNA
#'  with synthesis and degradation kinetics within the dynamic range of the experiment. In practice,
#'  you will (and should) be at least a bit lower than this as longer label times risk
#'  physiological impacts of metabolic labeling.
#'  \item Fraction labeled replicate correlation: This is the logit(fraction labeled)
#'  correlation between replicates, as inferred from your `metadf`.
#' }
#'
#' @param obj EZbakRData or EZbakRFractions object.
#' @param ... Parameters passed to the class-specific method.
#' If you have provided an EZbakRFractions object, then these can be (all play the
#' same role as in `EstimateKinetics()`, that is they get passed to `EZget()` to find
#' the fractions table you are interested in. See `?EstimateKinetics()` for details.):
#' \itemize{
#'  \item `features`
#'  \item `populations`
#'  \item `fraction_design`
#'
#' }
#' If you have provided an EZbakRData object, then these can be (all same the same
#' purpose as in `EstimateFractions`, so see `?EstimateFractions()` for details):
#' \itemize{
#'  \item `mutrate_populations`
#'  \item `features`
#'  \item `filter_cols`
#'  \item `filter_condition`
#'  \item `remove_features`
#' }
#' @import data.table
#' @importFrom magrittr %>%
#' @return A list of `ggplot2` objects visualizing the various aspects of your data
#' assessed by `EZQC()`.
#' @export
EZQC <- function(obj,
                 ...){

  UseMethod("EZQC")

}

#' Run quality control checks
#'
#' @param obj EZbakRFractions object, which is an EZbakRData object on which
#' `EstimateFractions` has been run.
#' @param ... Additional arguments. Currently goes unused.
#' @param features Set of features analyzed in the fractions table you are
#' interested QCing. This gets passed to `EZget()` to help find this table.
#' @param populations Set of mutation types analyzed in the fractions table
#' you are interested in QCing. This gets passed to `EZget()` to help find this table.
#' @param fraction_design The fraction "design matrix" specified to get the
#' fractions table you are interested in QCing. This gets passed to `EZget()` to
#' help find this table.
#' @import data.table
#' @importFrom magrittr %>%
#' @return A list of `ggplot2` objects visualizing the various aspects of your data
#' assessed by `EZQC()`.
#' @export
EZQC.EZbakRFractions <- function(obj,
                                 features = NULL,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 ...){

  ### Get mutation rate populations that were analyzed
  fname <- EZget(obj, type = "fractions",
                 features = features,
                 populations = populations,
                 fraction_design = fraction_design,
                 returnNameOnly = TRUE)
  mutrate_populations <- obj[["metadata"]][["fractions"]][[fname]][["populations"]]


  ### Check raw mutation rates
  graw <- check_raw_mutation_rates(obj,
                                   mutrate_populations = mutrate_populations)


  ### Check labeled and unlabeled mutation rates
  glabel <- check_plabeled(obj)


  ### Check read count replicate correlation
  gread <- check_read_count_corr_ezbf(obj,
                                      features = features,
                                      populations = populations,
                                      fraction_design = fraction_design)


  ### Check fraction labeled distribution
  gfld <- check_fl_dist(obj,
                        features = features,
                        populations = populations,
                        fraction_design = fraction_design)


  ### Check fraction labeled replicate correlation
  gflc <- check_fl_corr(obj,
                        features = features,
                        populations = populations,
                        fraction_design = fraction_design)


  return(
    list(
      Raw_mutrates = graw,
      Inferred_mutrates = glabel,
      Readcount_corr = gread,
      Fraction_labeled_dist = gfld,
      Fraction_labeled_corr = gflc
    )
  )

}



#' Run quality control checks
#'
#' @param obj EZbakRFractions object, which is an EZbakRData object on which
#' `EstimateFractions` has been run.
#' @param ... Additional arguments. Currently goes unused.
#' @param features Set of features analyzed in the fractions table you are
#' interested QCing. This gets passed to `EZget()` to help find this table.
#' @param populations Set of mutation types analyzed in the fractions table
#' you are interested in QCing. This gets passed to `EZget()` to help find this table.
#' @param fraction_design The fraction "design matrix" specified to get the
#' fractions table you are interested in QCing. This gets passed to `EZget()` to
#' help find this table.
#' @import data.table
#' @importFrom magrittr %>%
#' @return A list of `ggplot2` objects visualizing the various aspects of your data
#' assessed by `EZQC()`.
#' @export
EZQC.EZbakRArrowFractions <- function(obj,
                                 features = NULL,
                                 populations = NULL,
                                 fraction_design = NULL,
                                 ...){

  ### Get mutation rate populations that were analyzed
  fname <- EZget(obj, type = "fractions",
                 features = features,
                 populations = populations,
                 fraction_design = fraction_design,
                 returnNameOnly = TRUE)
  mutrate_populations <- obj[["metadata"]][["fractions"]][[fname]][["populations"]]


  ### Check raw mutation rates
  graw <- check_raw_mutation_rates(obj,
                                   mutrate_populations = mutrate_populations,
                                   arrow = TRUE)


  ### Check labeled and unlabeled mutation rates
  glabel <- check_plabeled(obj)


  ### Check read count replicate correlation
  gread <- check_read_count_corr_ezbf(obj,
                                      features = features,
                                      populations = populations,
                                      fraction_design = fraction_design)


  ### Check fraction labeled distribution
  gfld <- check_fl_dist(obj,
                        features = features,
                        populations = populations,
                        fraction_design = fraction_design)


  ### Check fraction labeled replicate correlation
  gflc <- check_fl_corr(obj,
                        features = features,
                        populations = populations,
                        fraction_design = fraction_design)


  return(
    list(
      Raw_mutrates = graw,
      Inferred_mutrates = glabel,
      Readcount_corr = gread,
      Fraction_labeled_dist = gfld,
      Fraction_labeled_corr = gflc
    )
  )

}



#' Run quality control checks
#'
#' @param obj An EZbakRData object
#' @param ... Additional arguments. Currently goes unused.
#' @param mutrate_populations Same as in `EstimateFractions()`. See `?EstimateFractions()`
#' for details.
#' @param features Same as in `EstimateFractions()`. See `?EstimateFractions()`
#' for details.
#' @param filter_cols Same as in `EstimateFractions()`. See `?EstimateFractions()`
#' for details.
#' @param filter_condition Same as in `EstimateFractions()`. See `?EstimateFractions()`
#' for details.
#' @param remove_features Same as in `EstimateFractions()`. See `?EstimateFractions()`
#' for details.
#' @import data.table
#' @import arrow
#' @importFrom magrittr %>%
#' @return A list of `ggplot2` objects visualizing the various aspects of your data
#' assessed by `EZQC()`.
#' @export
EZQC.EZbakRData <- function(obj,
                            mutrate_populations = "all",
                            features = "all",
                            filter_cols = "all",
                            filter_condition = `&`,
                            remove_features = c("NA", "__no_feature"),
                            ...){

  ### Check raw mutation rates
  graw <- check_raw_mutation_rates(obj,
                                   mutrate_populations = mutrate_populations)


  ### Check labeled and unlabeled mutation rates
  mutrates <- EstimateMutRates(obj,
                               populations = mutrate_populations)$mutation_rates


  glabel <- check_plabeled(obj,
                           mutrates)
  if(length(glabel) == 1){
    glabel <- glabel[[1]]
  }


  ### Check read count replicate correlation
  gread <- check_read_count_corr_ezbd(obj,
                                      features = features,
                                      filter_cols = filter_cols,
                                      filter_condition = filter_condition,
                                      remove_features = remove_features)



  return(
    list(
      Raw_mutrates = graw,
      Inferred_mutrates = glabel,
      Readcount_corr = gread
    )
  )

}


# Currently assumes cB, not arrow_dataset
check_raw_mutation_rates <- function(obj,
                                     mutrate_populations,
                                     arrow = FALSE){

  # Hack to address NOTEs
  n <- muttype <- mutrate <- NULL

  ### What mutation rates should be assessed?

  if(arrow){

    # cB columns
    cB <- obj$cBds
    cBschema <- arrow::schema(cB)
    cB_cols <- names(cBschema)

    ### Vectors of potential column names

    # Mutation counts possible
    mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                             c("T", "C", "G", "A", "U", "N"))
    mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

    illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

    mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]

    # Which mutcounts are in the cB?
    mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]


    # Base count columns in cB
    basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))


    # Which populations to analyze?
    if(mutrate_populations == "all"){

      muts_analyze <- mutcounts_in_cB

    }else{

      muts_analyze <- mutrate_populations

      if(!(muts_analyze %in% mutcounts_in_cB)){
        stop("You specified mutrate_populations not present in your cB!")
      }

    }

    nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


    ### Compute raw mutation rates of each type
    cB <- obj$cBds

    all_samples <- obj$metadf$sample

    mutrates <- vector(
      mode = "list",
      length = length(all_samples)
    )

    for(s in seq_along(all_samples)){


      cBsamp <- cBds %>%
        dplyr::collect()

      setDT(cBsamp)


      mutrates[[s]] <- cBsamp[, {
        rates <- sapply(seq_along(muts_analyze), function(i) {
          numerator <- sum(get(muts_analyze[i]) * n)
          denominator <- sum(get(nucs_analyze[i]) * n)
          rate <- numerator / denominator
          return(rate)
        })
        names(rates) <- paste0(muts_analyze)
        as.list(rates)
      }, by = sample]


    }

    mutrates <- dplyr::bind_rows(mutrates)


  }else{

    # Figure out which mutation counts are in the cB
    mutcounts_in_cB <- find_mutcounts(obj)

    basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

    # Which populations to analyze?
    if(mutrate_populations == "all"){

      muts_analyze <- mutcounts_in_cB

    }else{

      muts_analyze <- mutrate_populations

      if(!(muts_analyze %in% mutcounts_in_cB)){
        stop("You specified mutrate_populations not present in your cB!")
      }

    }

    nucs_analyze <- basecounts_in_cB[which(mutcounts_in_cB %in% muts_analyze)]


    ### Compute raw mutation rates of each type
    cB <- obj$cB

    mutrates <- cB[, {
      rates <- sapply(seq_along(muts_analyze), function(i) {
        numerator <- sum(get(muts_analyze[i]) * n)
        denominator <- sum(get(nucs_analyze[i]) * n)
        rate <- numerator / denominator
        return(rate)
      })
      names(rates) <- paste0(muts_analyze)
      as.list(rates)
    }, by = sample]


  }



  ### Make a plot of raw mutation rates of each type across samples
  glist <- vector(mode = "list",
                  length = length(muts_analyze))
  names(glist) <- muts_analyze

  # Pivot to make things easier when filtering and plotting in next step
  mutrates <- mutrates %>%
    tidyr::pivot_longer(
      cols = !!muts_analyze,
      names_to = "muttype",
      values_to = "mutrate"
    )

  for(m in seq_along(muts_analyze)){

    glist[[m]] <- mutrates %>%
      dplyr::filter(muttype == muts_analyze[m]) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = sample, y = mutrate)
      ) +
      ggplot2::geom_bar(stat = "identity",
                        position = "dodge",
                        fill = 'deepskyblue4') +
      ggplot2::theme_classic() +
      ggplot2::xlab("sample") +
      ggplot2::ylab("mutation rate") +
      ggplot2::labs(title = "Raw mutation rates") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90),
      )


  }

  if(length(glist) == 1){
    glist <- glist[[1]]
  }

  return(glist)

}




check_plabeled <- function(obj,
                           mutrates = NULL){


  # Hack to deal with NOTEs
  mutrate <- NULL


  # Gotta fetch mutrates
  if(is.null(mutrates)){

    mutrates <- obj$mutation_rates

  }

  mutrates_to_consider <- names(mutrates)[!grepl("_", names(mutrates))]

  # Plot distribution of pold and pnew for each sample
    # Filter out any potential hierarchical mutation rate estimates
  glist <- vector(mode = "list",
                  length = length(mutrates_to_consider))
  names(glist) <- mutrates_to_consider

  for(m in seq_along(mutrates_to_consider)){

    mutrate_subset <- mutrates[[mutrates_to_consider[m]]]

    glist[[mutrates_to_consider[m]]] <- mutrate_subset %>%
      tidyr::pivot_longer(
        cols = c("pold", "pnew"),
        names_to = "type",
        values_to = "mutrate"
      ) %>%
      dplyr::mutate(
        type = factor(type,
                      levels = c("pold", "pnew"))
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = sample, y = mutrate,
                     fill = type)
      ) +
      ggplot2::geom_bar(stat = "identity",
                        position = "dodge") +
      ggplot2::theme_classic() +
      ggplot2::xlab("sample") +
      ggplot2::ylab("mutation rate") +
      ggplot2::labs(title = "Labeled (red) and unlabeled (gray) read mutation rates",
                    subtitle = "Red bars ideally above blue line; Gray bars ideally below black line") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90),
        plot.title.position = "plot"
        ) +
      ggplot2::scale_fill_manual(values = c("darkgray", "darkred")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0.01,
                          linetype = "dotted",
                          size = 1.5,
                          color = "deepskyblue4") +
      ggplot2::geom_hline(yintercept = 0.004,
                          linetype = "dotted",
                          size = 1.5,
                          color = "black")


  }


  return(glist)

}



# Using EZbakRFractions object
check_read_count_corr_ezbf <- function(obj,
                                       features, populations,
                                       fraction_design){


  # Hack to deal with NOTEs
  `.` <- list()
  n <- NULL


  metadf <- data.table::copy(obj$metadf)

  rep_meta <- infer_replicates(obj,
                               consider_tl = FALSE)

  ### Filter out groups with no more than 1 replicate
  IDcol <- ifelse("replicate_id_imputed" %in% colnames(rep_meta),
                  "replicate_id_imputed",
                  "replicate_id")

  rep_meta <- rep_meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(IDcol))) %>%
    dplyr::filter(dplyr::n() > 1)

  if(nrow(rep_meta) == 0){
    message("Detected no conditions with replicates!
            Skipping read count replicate correlation analysis.")

    return(list())

  }


  ### Get read counts for each feature
  fname <- EZget(obj, type = "fractions",
                 features = features,
                 populations = populations,
                 fraction_design = fraction_design,
                 returnNameOnly = TRUE)
  fraction <- obj[["fractions"]][[fname]]
  features <- obj[["metadata"]][["fractions"]][[fname]]$features

  readcnts <- fraction %>%
    dplyr::filter(
      sample %in% rep_meta$sample
    ) %>%
    dplyr::select(sample, !!features, n) %>%
    dplyr::rename(reads = n)


  ### Make correlation plots

  glist <- make_corr_plots(readcnts,
                           rep_meta,
                           IDcol)


  return(glist)

}


# Using EZbakRData object
# Currently not compatible with arrow dataset
check_read_count_corr_ezbd <- function(obj,
                                       features,
                                       filter_cols,
                                       filter_condition,
                                       remove_features){


  # Hack to deal with NOTEs
  n <- NULL
  `.` <- list()


  cB <- data.table::copy(obj$cB)
  metadf <- data.table::copy(obj$metadf)

  rep_meta <- infer_replicates(obj,
                               consider_tl = FALSE)


  ### Filter out groups with no more than 1 replicate
  IDcol <- ifelse("replicate_id_imputed" %in% colnames(rep_meta),
                  "replicate_id_imputed",
                  "replicate_id")

  rep_meta <- rep_meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(IDcol))) %>%
    dplyr::filter(dplyr::n() > 1)

  if(nrow(rep_meta) == 0){
    message("Detected no conditions with replicates!
            Skipping read count replicate correlation analysis.")

    return(list())

  }


  ### Get read counts for each feature

  # What features are in cB?
  mutcounts_in_cB <- find_mutcounts(obj)

  basecounts_in_cB <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  cB_cols <- colnames(cB)

  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB,
                                             basecounts_in_cB,
                                             "sample", "n"))]

  # Need to determine which columns of the cB to group reads by
  if(features[1] == "all" & length(features) == 1){

    features_to_analyze <- features_in_cB

  }else{

    if(!all(features %in% features_in_cB)){

      stop("features includes columns that do not exist in your cB!")

    }else{

      features_to_analyze <- features

    }

  }

  group_cols <- c("sample", features_to_analyze)

  ### Filter out columns not mapping to any feature (easier and faster in data.table)

  if(filter_cols[1] == "all" & length(filter_cols) == 1){

    filter_cols <- features_to_analyze

  }

  cB <- cB[cB[, !Reduce(filter_condition, lapply(.SD, `%in%`, remove_features)),.SDcols = filter_cols], ]


  ### Get read counts

  readcnts <- cB[sample %in% rep_meta$sample][,.(reads = sum(n)),
                                             by = group_cols]


  ### Make correlation plots

  glist <- make_corr_plots(readcnts,
                           rep_meta,
                           IDcol)


  return(glist)

}



### MAKE READ COUNT CORRELATION PLOT
make_corr_plots <- function(table,
                                      rep_meta,
                                      IDcol,
                                      value = c("reads",
                                                "fraction"),
                                      fraction_type = "fraction_highTC",
                                      cutoff = 0.9){



  # Hack to deal with NOTEs
  density <- NULL


  value <- match.arg(value)

  axis_label_string <- ifelse(value == "reads",
                              "log10(read counts)",
                              paste0("logit(", fraction_type, ")"))


  RIDs <- unique(rep_meta[[IDcol]])

  glist <- vector(mode = "list",
                  length = length(RIDs))

  stat_df <- dplyr::tibble()
  for(r in seq_along(RIDs)){

    sub_meta <- rep_meta %>%
      dplyr::filter(!!dplyr::sym(IDcol) == RIDs[r])


    sub_table <- table %>%
      dplyr::filter(sample %in% sub_meta$sample)


    nreps <- nrow(sub_meta)

    subglist <- vector(mode = "list",
                       length = nreps)

    subplot_names <- c()
    count <- 1

    for(j in 1:(nreps - 1)){

      for(k in (j+1):nreps){


        samps <- c(sub_meta$sample[j],
                   sub_meta$sample[k])

        if(value == "reads"){

          corr_df <- sub_table %>%
            dplyr::filter(sample %in% samps) %>%
            dplyr::mutate(
              !!value := log10(!!dplyr::sym(value) + 1)
            ) %>%
            tidyr::pivot_wider(
              names_from = sample,
              values_from = !!value
            )

        }else{

          corr_df <- sub_table %>%
            dplyr::filter(sample %in% samps) %>%
            dplyr::mutate(
              !!value := logit(!!dplyr::sym(value))
            ) %>%
            tidyr::pivot_wider(
              names_from = sample,
              values_from = !!value
            )

        }


        subglist[[count]] <- corr_df %>%
          stats::na.omit() %>%
          dplyr::mutate(
            density = get_density(
              x = !!dplyr::sym(samps[1]),
              y = !!dplyr::sym(samps[2]),
              n = 200
            )
          ) %>%
          stats::na.omit() %>%
          ggplot2::ggplot(
            aes(x = .data[[samps[1]]],
                y = .data[[samps[2]]],
                color = density)
          ) +
          ggplot2::geom_point() +
          ggplot2::theme_classic() +
          ggplot2::scale_color_viridis_c() +
          ggplot2::geom_abline(
            slope = 1,
            intercept = 0,
            color = 'darkred',
            linetype = 'dotted'
          ) +
          ggplot2::xlab(paste0(samps[1], " ", axis_label_string)) +
          ggplot2::ylab(paste0(samps[2], " ", axis_label_string))


        stat_df <- stat_df %>%
          dplyr::bind_rows(
            dplyr::tibble(
              sample_1 = samps[1],
              sample_2 = samps[2],
              correlation = stats::cor(corr_df[[samps[1]]],
                                corr_df[[samps[2]]], use = "complete.obs")
            )
          )

        subplot_names[count] <- paste0(samps[1], "_vs_", samps[2])
        count <- count + 1

      }



    }

    names(subglist) <- subplot_names
    glist[[r]] <- subglist

  }

  ### Print correlations
  message(paste0(c(paste0(axis_label_string, " correlation for each pair of replicates are:"), utils::capture.output(stat_df)), collapse = "\n"))


  message("")
  if(any(stat_df$correlation < cutoff)){
    message(paste0(axis_label_string, " correlation is low in one or more samples.
              Did you properly specify all sample detail columns in your metadf?")
            )
  }else{
    message(paste0(axis_label_string, " correlations are high, suggesting good reproducibility!"))
  }
  message("")

  names(glist) <- paste0("Group_", 1:length(glist))

  return(glist)

}





check_fl_dist <- function(obj,
                          features,
                          populations,
                          fraction_design){

  # Hack to address NOTEs
  fraction <- NULL

  ### Get read counts for each feature
  fname <- EZget(obj, type = "fractions",
                 returnNameOnly = TRUE)
  fraction_table <- obj[["fractions"]][[fname]]
  features <- obj[["metadata"]][["fractions"]][[fname]]$features

  fraction_cols <- colnames(fraction_table)
  fraction_cols <- fraction_cols[grepl("^fraction_", fraction_cols)]


  glist <- vector(mode = "list",
                  length = length(fraction_cols))

  avg_fractions <- dplyr::tibble()

  for(fc in seq_along(fraction_cols)){

    fractions <- fraction_table %>%
      dplyr::select(sample, !!features, !!fraction_cols[fc]) %>%
      dplyr::rename(fraction = !!dplyr::sym(fraction_cols[fc])) %>%
      dplyr::filter(fraction != 0)

    avg_fractions <- avg_fractions %>%
      dplyr::bind_rows(
        fractions %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(
            avg_fraction = mean(fraction)
          ) %>%
          dplyr::mutate(
            fraction_type = fraction_cols[fc]
          )
      )


    samples <- unique(fractions$sample)

    for(s in seq_along(samples)){



      glist[[fc]][[samples[s]]] <- fractions %>%
        ggplot2::ggplot(aes(x = fraction)) +
        ggplot2::geom_density() +
        ggplot2::theme_classic() +
        ggplot2::xlab(fraction_cols[fc]) +
        ggplot2::ylab("Density")

    }

  }

  message(paste0(c("Average fractions for each sample are:", utils::capture.output(avg_fractions)), collapse = "\n"))


  lower_fxn_cutoff <- ifelse(
    length(fraction_cols) == 1,
    0.5/10,
    1/(10*length(fraction_cols))
  )
  upper_fxn_cutoff <- 1 - lower_fxn_cutoff


  message("")
  low_fx <- FALSE
  high_fx <- FALSE
  if(any(avg_fractions$avg_fraction < lower_fxn_cutoff)){
    message("One or more of your samples have very low amounts of labeling.
            This could mean that your label times are short, which can limit
            the statistical power of EZbakR analyses.")
    low_fx <- TRUE
  }

  if(any(avg_fractions$avg_fraction > upper_fxn_cutoff)){
    message("One or more of your samples have very high amounts of labeling.
            This could mean that your label times are long, which can limit
            the statistical power of EZbakR analyses.")
    high_fx <- TRUE
  }

  if(!low_fx & !high_fx){
    message("Labeling rates (e.g., fraction labeled for single label experiments) look good!")
  }

  message("")

  if(length(glist) == 1){

    glist <- glist[[1]]

  }else{

    names(glist) <- fraction_cols

  }


  return(glist)

}

check_fl_corr <- function(obj,
                          features,
                          populations,
                          fraction_design){

  metadf <- data.table::copy(obj$metadf)

  rep_meta <- infer_replicates(obj,
                               consider_tl = TRUE)

  ### Filter out groups with no more than 1 replicate
  IDcol <- ifelse("replicate_id_imputed" %in% colnames(rep_meta),
                  "replicate_id_imputed",
                  "replicate_id")

  rep_meta <- rep_meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(IDcol))) %>%
    dplyr::filter(dplyr::n() > 1)

  if(nrow(rep_meta) == 0){
    message("Detected no conditions with replicates!
            Skipping read count replicate correlation analysis.")

    return(list())

  }


  ### Get read counts for each feature
  fname <- EZget(obj, type = "fractions",
                 returnNameOnly = TRUE)
  fraction_table <- obj[["fractions"]][[fname]]
  features <- obj[["metadata"]][["fractions"]][[fname]]$features

  fraction_cols <- colnames(fraction_table)
  fraction_cols <- fraction_cols[grepl("^fraction_", fraction_cols)]


  glist <- vector(mode = "list",
                  length = length(fraction_cols))

  for(fc in seq_along(fraction_cols)){

    fractions <- fraction_table %>%
      dplyr::filter(
        sample %in% rep_meta$sample
      ) %>%
      dplyr::select(sample, !!features, !!fraction_cols[fc]) %>%
      dplyr::rename(fraction = !!dplyr::sym(fraction_cols[fc]))


    glist[[fc]] <- make_corr_plots(
      setDT(data.table::copy(fractions)),
      rep_meta,
      IDcol,
      value = "fraction",
      fraction_type = fraction_cols[fc],
    )

  }

  if(length(glist) == 1){

    glist <- glist[[1]]

  }else{

    names(glist) <- fraction_cols

  }


  return(glist)
}


# Replicate = group of samples with same sample details
# in metadf.
infer_replicates <- function(obj,
                             consider_tl){

  metadf <- data.table::copy(obj$metadf)

  if("replicate_id" %in% colnames(metadf)){
    new_col <- "replicate_id_imputed"
  }else{
    new_col <- "replicate_id"
  }

  # Find label time columns
  mutcounts_in_cB <- find_mutcounts(obj)

  tl_cols_possible <- c("tl", "tpulse", "tchase",
                        paste0("tl_", mutcounts_in_cB),
                        paste0("tpulse_", mutcounts_in_cB),
                        paste0("tchase_", mutcounts_in_cB))


  tl_cols <- tl_cols_possible[tl_cols_possible %in% colnames(metadf)]


  if(!consider_tl){

    cols_to_filter <- c("sample", tl_cols)

  }else{

    # Need to remove label free samples
    metadf <- metadf %>%
      dplyr::rowwise() %>%
      dplyr::filter(!all(dplyr::c_across(dplyr::all_of(tl_cols)) == 0))


    cols_to_filter <- "sample"

  }

  # Filter and add a replicate ID column
  replicates <- metadf %>%
    dplyr::ungroup() %>%
    dplyr::select(-!!cols_to_filter) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!new_col := 1:dplyr::n())



  rep_meta <- metadf %>%
    dplyr::inner_join(replicates,
                      by = colnames(replicates)[colnames(replicates) != new_col])

  return(rep_meta)

}



### Get point density so that I can color by density
### Source: https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



### Some helper functions borrowed from MASS
### Due to MASS now requiring R > 4.4, I'd like to avoid
### a strict requirement of MASS. Technically ggplot2, which I
### import in EZbakR, requires MASS, but this might change in the
### near future (https://github.com/tidyverse/ggplot2/issues/5986).
### I could also suggest MASS, check if it is installed, and then
### avoid coloring by density if not...
### Source for code: https://github.com/cran/MASS/tree/master/R
bandwidth.nrd <- function (x)
{
  r <- stats::quantile(x, c(0.25, 0.75))
  h <- (r[2L] - r[1L])/1.34
  4 * 1.06 * min(sqrt(stats::var(x)), h) * length(x)^(-1/5)
}

kde2d <- function (x, y, h, n = 25, lims = c(range(x), range(y)))
{
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  if (any(!is.finite(x)) || any(!is.finite(y)))
    stop("missing or infinite values in the data are not allowed")
  if (any(!is.finite(lims)))
    stop("only finite values are allowed in 'lims'")
  n <- rep(n, length.out = 2L)
  gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
  h <- if (missing(h))
    c(bandwidth.nrd(x), bandwidth.nrd(y))
  else rep(h, length.out = 2L)
  if (any(h <= 0))
    stop("bandwidths must be strictly positive")
  h <- h/4
  ax <- outer(gx, x, "-")/h[1L]
  ay <- outer(gy, y, "-")/h[2L]
  z <- tcrossprod(matrix(stats::dnorm(ax), , nx), matrix(stats::dnorm(ay),
                                                  , nx))/(nx * h[1L] * h[2L])
  list(x = gx, y = gy, z = z)
}
