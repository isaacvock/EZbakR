#' `EZbakRDataobject` constructor for internal use
#'
#' \code{new_EZbakRData} efficiently creates an object of class `EZbakRData`.
#' It does not perform any rigorous checks of the legitimacy of this object.
#' @param cB Data frame tracking the sample ID, mutational and nucleotide content,
#' and feature assignment of sequencing reads.
#' @param metadf Data frame tracking features of each of the samples included
#' in `cB`.
new_EZbakRData <- function(cB, metadf){
  stopifnot(is.data.frame(cB))
  stopifnot(is.data.frame(metadf))
  structure(list(cB = cB, metadf = metadf), class = "EZbakRData")
}

#' `EZbakRDataobject` validator
#'
#' \code{validate_EZbakRData} ensures that input for `EZbakRData` construction
#' is valid.
#' @import data.table
#' @param obj An object of class `EZbakRData`
validate_EZbakRData <- function(obj){

  ### Vectors of potential column names

  # Mutation counts
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[,1], mutcounts[,2])

  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")

  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]


  ### Extract objects and info from objects

  vals <- unclass(obj)

  cB <- dplyr::as_tibble(vals$cB)
  metadf <- dplyr::as_tibble(vals$metadf)

  cB_cols <- colnames(cB)
  metadf_cols <- colnames(metadf)

  # Which mutcounts are in the cB?
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]


  ### Does cB_cols contain sample and n?
  if(!("sample" %in% cB_cols & "n" %in% cB_cols)){

    rlang::abort(
      "cB must include columns named sample (representing the sample ID)
      and n (representing the number of reads with identical values for the
      columns included in the cB)!",
      class = "cB_sample_n"
    )

  }


  ### Does cB_cols contain the necessary base count column(s)?
  basecounts_expected <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))

  if(!all(basecounts_expected %in% cB_cols)){

    rlang::abort(
      "Not all of the columns tracking the counts of relevant nucleotides
      are present! For example, if you have columns called TC and GA,
      you need to have corresponding columns called nT and nG. TC and GA
      should represent the number of T-to-C and G-to-A mutations
      in sequencing reads, respectively, and nT and nG should represent
      the number of Ts and Gs in the read (if there had been no mutations in that read).",
      class = "cB_basecount"
    )

  }


  ### Does metadf contain sample and correct tl(s)?

  if(!("sample" %in% metadf_cols)){
    rlang::abort(
      "metadf must include a column named sample (representing the sample ID)!",
      class = "metadf_sample"
    )
  }

  # If only one mutation type tracked in cB, then need either tl or tchase + tpulse
  # If more than one mutation type tracked in cB, need tl_<muttype> for each
  # <muttype (e.g., TC, GA, etc.) in the cB.
  if(length(mutcounts_in_cB) > 1){

    tl_expected <- paste0("tl_", mutcounts_in_cB)
    tl_expected_p <- paste0("tpulse_", mutcounts_in_cB)
    tl_expected_c <- paste0("tchase_", mutcounts_in_cB)

    if(!all(tl_expected %in% metadf_cols)){

      if(!(all(tl_expected_p %in% metadf_cols) & all(tl_expected_c %in% metadf_cols))){


        rlang::abort(
          "Not all of the relevant label times are included in the metadf.
          For example, if your cB has columns TC and GA, your metadf must
          include columns called tl_TC and tl_GA, representing the label times
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




  ### Check if all of the samples in cB are also in metadf

  samps_cB <- unique(cB$sample)

  ### TO-DO: COULD THROW A WARNING IF THE INVERSE IS TRUE.
  if(!all(metadf$sample %in% samps_cB)){
    rlang::abort("Not all samples in the metadf are present in the cB!",
          class = "cB_metadf_samples")
  }


  ### Check that mutation counts can be coerced to positive integers
  is_pos_whole <- function(x){

    if(all(is.numeric(x))){

      bool <- all((floor(x) == x) & (x >= 0))

    }else{

      bool <- FALSE

    }

    return(bool)


  }
  is_pos_whole2 <- function(x){

    if(all(is.numeric(x))){

      bool <- all((floor(x) == x) & (x > 0))

    }else{

      bool <- FALSE

    }

    return(bool)

  }

  if(!all(sapply(cB[,mutcounts_in_cB], is_pos_whole))){

    rlang::abort("Not all columns tracking counts of mutations are positive whole numbers!",
          class = "mut_not_pos_whole")

  }

  ### Check that base counts can be coerced to positive integers

  if(!all(sapply(cB[,basecounts_expected], is_pos_whole))){

    rlang::abort("Not all columns tracking counts of nucleotides are positive whole numbers!",
                 class = "basecounts_not_pos_whole")

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


  ### Check that n in cB is numeric and > 0

  if(!all(sapply(cB[,"n"], is_pos_whole2))){

    rlang::abort("Not all columns of cB tracking counts of reads contain positive whole numbers
         strictly greating than 0!",
          class = "cB_n_pos_whole")

  }


  ### Make sure there is at least one labeled sample
    ## Also label for each mutation type tracking; tell them to get rid
    ## of unnecessary mutation types if not.
    ## Also make sure for pulse-chase that there is in fact a pulse and a
    ## chase


  ### Make sure there is at least one feature column in cB


  ### Convert all other metadf columns to factors

  cols_to_convert <- metadf_cols[!(metadf_cols %in% c(tl_cols, "sample"))]

  metadf <- as.data.frame(metadf)

  metadf[cols_to_convert] <- lapply(metadf[cols_to_convert], as.factor)

  metadf <- dplyr::as_tibble(metadf)

  obj$metadf <- metadf

  return(obj)




}


#' `EZbakRData` object helper function for users
#'
#' \code{EZbakRData} creates an object of class `EZbakRData` and checks the validity
#' of the provided input.
#' @param cB Data frame with the following columns:
#' \itemize{
#'  \item sample: Name given to particular sample from which data was collected.
#'  \item mutational counts: Integers corresponding to the number of a particular
#'  mutation seen in a sequencing read. The following column names are allowed:
#'  \itemize{
#'    \item TC: Number of Thymine-to-Cytosine mutations
#'    \item TA: Number of Thymine-to-Adenine mutations
#'    \item TG: Number of Thymine-to-Guanine mutations
#'    \item CT: Number of Cytosine-to-Thymine mutations
#'    \item CA: Number of Cytosine-to-Adenine mutations
#'    \item CG: Number of Cytosine-to-Guanine mutations
#'    \item CU: Number of Cytosine-to-Uridine mutations
#'    \item AT: Number of Adenine-to-Thymine mutations
#'    \item AC: Number of Adenine-to-Cytosine mutations
#'    \item AG: Number of Adenine-to-Guanine mutations
#'    \item AU: Number of Adenine-to-Uridine mutations
#'    \item GT: Number of Guanine-to-Thymine mutations
#'    \item GC: Number of Guanine-to-Cytosine mutations
#'    \item GA: Number of Guanine-to-Adenine mutations
#'    \item GU: Number of Guanine-to-Uridine mutations
#'    \item TN: Number of Thymine-to-Adenine/Cytosine/Guanine mutations
#'    \item CN: Number of Cytosine-to-Adenine/Thymine/Guanine/Uridine mutations
#'    \item AN: Number of Adenine-to-Thymine/Cytosine/Guanine/Uridine mutations
#'    \item GN: Number of Guanine-to-Adenine/Cytosine/Thymine/Uridine mutations
#'    \item UN: Number of Uridine-to-Adenine/Cytosine/Guanine mutations
#'    \item NT: Number of Adenine/Cytosine/Guanine-to-Thymine mutations
#'    \item NC: Number of Adenine/Thymine/Guanine/Uridine-to-Cytosine mutations
#'    \item NtoA: Number of Thymine/Cytosine/Guanine/Uridine-to-Adenine mutations. (Naming convention changed because NA taken)
#'    \item NU: Number of Cytosine/Guanine/Adenine-to-Uridine mutations.
#'    \item NN: Number of any kind of mutation
#'  }
#'  \item base nucleotide count: Integers corresponding to the number of instances of a
#'  particular type of nucleotide whose mutations are tracked in a corresponding
#'  mutation count column. The following column names are allowed:
#'  \itemize{
#'    \item nT: Number of Thymines
#'    \item nG: Number of Guanines
#'    \item nA: Number of Adenines
#'    \item nC: Number of Cytosines
#'    \item nU: Number of Uridines
#'    \item nN: number of any kind of nucleotide
#'  }
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
#'
#' }
#' @return An EZbakRData object. This is simply a list of the provide `cB` and
#' `metadf` with class `EZbakRData`
#' @examples
#'
#' # Simulate data
#' simdata <- EZSimulate(30)
#'
#' # Create EZbakRData object
#' ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
#'
#' @import data.table
#' @export
EZbakRData <- function(cB, metadf){

  cB <- setDT(cB)
  metadf <- setDT(metadf)

  validate_EZbakRData(new_EZbakRData(cB, metadf))

}
