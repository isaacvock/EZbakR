#' `EZbakRFractions` object constructor
#'
#' \code{new_EZbakRData} efficiently creates an object of class `EZbakRData`.
#' It does not perform any rigorous checks of the legitimacy of this object.
#' @param fractions Data frame containing information about the fraction of reads
#' from each mutational population of interest.
#' @param metadf Data frame tracking features of each of the samples included
#' in `cB`.
new_EZbakRFractions <- function(fractions, metadf){
  stopifnot(is.data.frame(fractions))
  stopifnot(is.data.frame(metadf))
  structure(list(cB = fractions, metadf = metadf), class = "EZbakRFractions")
}

#' `EZbakRFractions` object validator
#'
#' \code{validate_EZbakRData} ensures that input for `EZbakRFractions` construction
#' is valid.
#' @param obj An object of class `EZbakRFractions`
validate_EZbakRFractions <- function(obj){

  ### Vector of potential mutational populations

  # Mutation counts
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]


  ### Extract objects and info from objects

  vals <- unclass(obj)

  fractions <- dplyr::as_tibble(vals$fractions)
  metadf <- dplyr::as_tibble(vals$metadf)

  fraction_cols <- colnames(fractions)
  metadf_cols <- colnames(metadf)


  ### Get the mutation counts of interest

  # A bit sketchy because you risk snagging a feature column,
  # but it would have to specifcally be a feature column that includes 'fraction_'
  # in its name, so that seems like an unlikely occurence. Furthermore,
  # a substring of the feature column would have to include 'high' or 'low'
  # to screw things up. Even then, I will probably add a filter for valid
  # mutation types to catch this crazy edge case (or a misnamed fraction).
  substrings <- unlist(strsplit(fraction_cols[grepl("fraction_", fraction_cols)], "_"))
  pops <- substrings[grepl("high", substrings) | grepl("low", fraction_cols)]
  mutations_in_fractions <- substr(pops, nchar(pops) - 1, nchar(pops))

  mutations_in_fractions <- mutations_in_fractions[mutations_in_fractions %in% mutcounts]


  ### Does fractions data frame contain columns named "sample" and "n"?
  if(!("sample" %in% fraction_cols & "n" %in% fraction_cols)){

    rlang::abort(
      "fractions must include columns named sample (representing the sample ID)
      and n (representing the number of reads)!",
      class = "fraction_sample_n"
    )

  }



  ### Does metadf contain sample and correct tl(s)?

  if(!("sample" %in% metadf_cols)){
    rlang::abort(
      "metadf must include a column named sample (representing the sample ID)!",
      class = "metadf_sample"
    )
  }

  # If only one mutation type tracked in fractions, then need either tl or tchase + tpulse
  # If more than one mutation type tracked in fractions, need tl_<muttype> for each
  # <muttype (e.g., TC, GA, etc.) in the fractions.
  if(length(mutations_in_fractions) > 1){

    tl_expected <- paste0("tl_", mutations_in_fractions)
    tl_expected_p <- paste0("tpulse_", mutations_in_fractions)
    tl_expected_c <- paste0("tchase_", mutations_in_fractions)

    if(!all(tl_expected %in% metadf_cols)){

      if(!(all(tl_expected_p %in% metadf_cols) & all(tl_expected_c %in% metadf_cols))){


        rlang::abort(
          "Not all of the relevant label times are included in the metadf.
          For example, if your fractions has columns referring to high/low TC and GA,
          your metadf must include columns called tl_TC and tl_GA, representing the label times
          for metabolic labels whose recoding yield apparent T-to-C and
          G-to-A mutations, respectively. Alternatively, if you performed a
          pulse-chase, then you need columns tpulse_TC, tpulse_GA, tchase_TC,
          and tpulse_GA.",
          class = "metadf_tl_multimut"
        )

      }else{

        tl_cols <- c(tl_expected_p, tl_expected_c)

      }



    }else{

      tl_cols <- tl_expected

    }




  }else{

    if(!("tl" %in% metadf_cols)){

      if(!all(c("tchase", "tpulse") %in% metadf_cols)){

        rlang::abort(
          "metadf must contain a column called tl (the label time if using a pulse-label design)
          or a combination of tpulse and tchase (pulse and chase times respectively
          if using a pulse-chase design).",
          class = "metadf_tl_onemut"
        )

      }else{

        tl_cols <- c("tchase", "tpulse")

      }

    }else{

      tl_cols <- "tl"

    }

  }




  ### Check if all of the samples in fractions are also in metadf

  samps_fractions <- unique(fractions$sample)

  if(!all(samps_fractions %in% metadf$sample)){
    rlang::abort("Not all samples in fractions are present in the metadf!",
                 class = "fractions_metadf_samples")
  }


  ### Check that n in fractions is numeric and > 0

  is_pos_whole2 <- function(x){

    if(all(is.numeric(x))){

      bool <- all((floor(x) == x) & (x > 0))

    }else{

      bool <- FALSE

    }

    return(bool)

  }


  if(!all(sapply(fractions[,"n"], is_pos_whole2))){

    rlang::abort("Not all columns of fraction tracking counts of reads contain positive whole numbers
         strictly greating than 0!",
                 class = "fractions_n_pos_whole")

  }


  ### Check that label times are numeric and >= 0


  is_pos_num <- function(x){

    if(all(is.numeric(x))){

      bool <- all(x >= 0)


    }else{

      bool <- FALSE

    }


    return(bool)

  }


  if(!all(sapply(metadf[,tl_cols], is_pos_num))){

    rlang::abort("Not all columns of metadf representing label times are numbers
         greater than or equal to 0.",
                 class = "tl_pos_num")

  }





  ### Make sure there is at least one labeled sample
  ## Also label for each mutation type tracking; tell them to get rid
  ## of unnecessary mutation types if not.
  ## Also make sure for pulse-chase that there is in fact a pulse and a
  ## chase


  ### Make sure there is at least one feature column in fractions


  return(obj)




}


#' `EZbakRFractions` helper function for users
#'
#' \code{EZbakRData} creates an object of class `EZbakRData` and checks the validity
#' of the provided input.
#' @param fractions Data frame with the following columns:
#' \itemize{
#'  \item sample: Name given to particular sample from which data was collected.
#'  \item estimates of population fractions: These columns refer to the estimate
#'  for the fraction of reads coming from a particular mutational population.
#'  For example, in a standard NR-seq experiment, you should have one column
#'  named `fraction_highTC`. This refers to the fraction of RNA inferred to have
#'  a high T-to-C mutation rate (e.g., the newly synthesized RNA in a pulse-labeling NR-seq
#'  experiment).
#'  If you estimated the fractions of more than 2 mutation types (e.g., T-to-C and
#'  G-to-A), then you need to explicitly list all fractions of interest estimated.
#'  For example, in a TILAC experiment, this would be `fraction_highTC_lowGA`,
#'  `fraction_lowTC_highGA`, and `fraction_lowTC_lowGA`.
#'  \item features: Any columns that cannot be interpreted as a mutation count
#'  or base nucleotide count (and that aren't named `sample` or `n`) will be
#'  interpreted as an ID for a genomic "feature" from which a read originated.
#'  Common examples of features and typical column names for said features include:
#'  \itemize{
#'    \item Genes; common column names: gene, gene_id, gene_name, GF
#'    \item Genes-exonic; common column names: gene_exon, gene_id_exon, gene_name_exon, XF
#'    \item Transcripts; common column names: transcripts, TF
#'    \item Exonic bins; common column names: exonic_bins, EF, EB
#'    \item Exons; common column names: exons, exon_ids
#'  }
#'  In some cases, a read will often map to multiple features (e.g., exons). Many
#'  functions in bakR expect each of the feature IDs in these cases to be separated
#'  by `+`. For example, if a read overlaps with two exons, with IDs exon_1 and exon_2,
#'  then the corresponding entry in a  column of exonic assignments would be "exon_1+exon_2".
#'  The default expectation can be overwritten though and is thus not strictly enforced.
#'  \item n: Number of reads with identical values for all other columns.
#' }
#' @param metadf Data frame detailing various aspects of each of the samples included
#' in the cB. This includes:
#' \itemize{
#'  \item `sample`: The sample ID, which should correspond to a sample ID in the provided cB.
#'  \item `tl`: Metabolic label time. There are several edge cases to be aware of:
#'  \itemize{
#'    \item If more than one metabolic label was used in the set of samples described
#'    by the metadf (e.g., s4U and s6G were used), then the `tl` column should be
#'    replaced by `tl_<muttype>`, where `<muttype>` represents the corresponding mutation
#'    type count column in the cB that the label whose incubation time will be listed
#'    in this column. For example, if feeding with s4U in some samples and s6G in others,
#'    then performing standard nucleotide recoding chemistry, you will include
#'    `tl_TC` and `tl_GA` columns corresponding to the s4U and s6G label times, respectively.
#'    \item If a pulse-chase experimental design was used (!!this is strongly discouraged
#'    unless you have a legitimate reason to prefer this design to a pulse-label
#'    design!!), then you should have columns named `tpulse` and `tchase`, corresponding
#'    to the pulse and chase times respectively. The same _<muttype> convention should
#'    be used in the case of multi-label pulse-chase designs.
#'  }
#' \item sample characteristics: The remaining columns can be named whatever you like
#' and should include distinguishing features of groups of samples. Common columns might
#' include:
#' \itemize{
#'  \item `treatment`: The experimental treatment applied to a set of samples.
#'  This could represent things like genetic knockouts or knockdowns, drug treatments, etc.
#'  \item `batch`: An ID for sets of samples that were collected and/or processed together.
#'  Useful for regressing out technical batch effects
#'  }
#' \item `assay`: This optional column should include a string that
#' describes the type of experiment that was done so as to influence
#' how EZbakR analyzes and interprets the data from those samples.
#' Possible values for `assay` currently include:
#' \itemize{
#'  \item standard: Refers to the "standard" nucleotide recoding RNA-seq methods
#'  (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.), in which cells are fed with a
#'  single metabolic label, RNA is extracted and sequenced, and mutations of a particular
#'  type are counted
#'  \item STL: Refers to Start-TimeLapse-seq, a method combining Start-seq (developed
#'  by Karen Adelman's lab) with TimeLapse-seq. Used to infer the kinetics of
#'  transcription initiation and promoter-proximal pause site departure.
#'  \item TT: Refers to Transient-Transcriptome NR-seq, a method combining TT-seq (developed
#'  by Patrick Cramer's lab) with NR-seq. TT-seq involves biochemically enriching
#'  for labeled RNA. By combining this method with nucleotide recoding chemistry (as was
#'  first done by the Simon lab with TT-TimeLapse-seq and has since been done with SLAM
#'  chemistry, often referred to as TTchem-seq), it is possible to bioinformatically
#'  filter out reads coming from unlabeled RNA background.
#'  \item TILAC: Refers to TILAC, a method developed by the Simon lab to achieve
#'  spike-in free normalization of RNA-seq data through the use of a dual labeling
#'  approach inspired by the proteomic method SILAC.
#'  \item subcellular: Refers to techniques such as subcellular TimeLapse-seq (developed
#'  by Stirling Churchman's lab) which combine subcellular fractionation with
#'  NR-seq to infer additional kinetic parameters.
#'  \item sc: Refers to single-cell RNA-seq implementations of NR-seq.
#'  }
#'
#' }
#' @return An EZbakRData object. This is simply a list of the provide `cB` and
#' `metadf` with class `EZbakRData`
#' @export
EZbakRFractions <- function(fractions, metadf){

  fractions <- setDT(fractions)
  metadf <- setDT(metadf)

  validate_EZbakRFractions(new_EZbakRFractions(fractions, metadf))

}