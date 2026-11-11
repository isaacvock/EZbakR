# Internal helper: canonical vector of valid mutation count column names
ezbakr_valid_mutcounts <- function(){
  mutcounts <- expand.grid(c("T", "C", "G", "A", "U", "N"),
                           c("T", "C", "G", "A", "U", "N"))
  mutcounts <- paste0(mutcounts[, 1], mutcounts[, 2])
  illegal_mutcounts <- c("TT", "CC", "GG", "AA", "UU")
  mutcounts <- mutcounts[!(mutcounts %in% illegal_mutcounts)]
  mutcounts[mutcounts == "NA"] <- "NtoA"
  mutcounts
}


# Internal helper: unwrap nested arrow_dplyr_query objects to their root Dataset
unwrap_arrow_dataset <- function(x){
  base_ds <- x
  while(inherits(base_ds, "arrow_dplyr_query")){
    base_ds <- base_ds$.data
  }
  base_ds
}


# Internal helper: strip R-specific schema metadata (e.g., legacy data.table
# `.internal.selfref` externalptr attributes in Parquet footers) via WithMetadata(NULL)
# without deserializing `schema$metadata`, preventing Arrow unsafe-metadata warnings.
strip_arrow_r_metadata <- function(x){
  if(is.null(x)) return(x)
  if(inherits(x, "arrow_dplyr_query")){
    x$.data <- strip_arrow_r_metadata(x$.data)
    return(x)
  }
  if(methods::is(x, "Dataset")){
    x <- tryCatch(
      x$WithSchema(x$schema$WithMetadata(NULL)),
      error = function(e) x
    )
    return(x)
  }
  x
}


# Internal helper: check whether an Arrow Dataset / query has live C++ pointers
arrow_is_live <- function(x){
  if(is.null(x)) return(FALSE)
  base_ds <- unwrap_arrow_dataset(x)
  if(is.null(base_ds) || !(methods::is(base_ds, "ArrowObject") && methods::is(base_ds, "Dataset"))){
    return(FALSE)
  }
  tryCatch({
    sc_base <- arrow::schema(base_ds)
    sc_top <- arrow::schema(x)
    # Also touch base_ds$files to verify the FileSystemDataset pointer itself is valid
    fls <- base_ds$files
    !is.null(sc_base) && !is.null(sc_top) && is.character(fls)
  }, error = function(e) FALSE)
}


# Internal helper: infer file format from file paths when not explicitly supplied
infer_arrow_format_from_paths <- function(paths){
  probe <- paths
  if(length(paths) == 1L && dir.exists(paths)){
    probe <- list.files(paths, recursive = TRUE, full.names = TRUE)
  }
  if(length(probe) == 0L) return("parquet")
  low <- tolower(basename(probe))
  if(any(grepl("\\.(parquet|pq)$", low))) return("parquet")
  if(any(grepl("\\.(feather|arrow|ipc)$", low))) return("feather")
  if(any(grepl("\\.tsv(\\.gz)?$", low))) return("tsv")
  if(any(grepl("\\.csv(\\.gz)?$", low))) return("csv")
  "parquet"
}


# Internal helper: resolve a path, file vector, Dataset, or arrow_dplyr_query
resolve_arrow_input <- function(cBds,
                                format = NULL,
                                partitioning = NULL,
                                hive_style = TRUE,
                                ...){
  if(is.character(cBds)){
    if(length(cBds) == 0L || anyNA(cBds) || any(!nzchar(cBds))){
      rlang::abort(
        "Path(s) supplied to `cBds` must be non-empty character strings!",
        class = "arrow_path_invalid"
      )
    }
    exists_mask <- file.exists(cBds) | dir.exists(cBds)
    if(!all(exists_mask)){
      missing_p <- cBds[!exists_mask]
      rlang::abort(
        paste0("Arrow dataset path(s) not found on disk: ",
               paste(sprintf("'%s'", missing_p), collapse = ", ")),
        class = "arrow_path_not_found"
      )
    }

    norm_paths <- normalizePath(cBds, winslash = "/", mustWork = TRUE)
    fmt <- if(is.null(format)) infer_arrow_format_from_paths(norm_paths) else format
    is_single_dir <- (length(norm_paths) == 1L) && dir.exists(norm_paths)
    dots <- list(...)

    if(is_single_dir){
      open_args <- list(sources = norm_paths, format = fmt, hive_style = hive_style)
      if(!is.null(partitioning)){
        open_args$partitioning <- partitioning
      }
      if(length(dots) > 0L){
        open_args <- c(open_args, dots)
      }
      ds <- do.call(arrow::open_dataset, open_args)
    }else{
      ds <- open_arrow_files_source(norm_paths,
                                    fmt = fmt,
                                    partitioning = partitioning,
                                    hive_style = hive_style,
                                    dots = dots)
    }
    ds <- strip_arrow_r_metadata(ds)
    return(list(
      ds = ds,
      source_type = if(is_single_dir) "directory" else "files",
      root_path = if(is_single_dir) norm_paths else NULL,
      explicit_paths = norm_paths,
      format = fmt,
      partitioning_cols = if(is.character(partitioning)) partitioning else NULL,
      hive_style = hive_style,
      open_dots = dots
    ))
  }

  if(methods::is(cBds, "Dataset") || inherits(cBds, "arrow_dplyr_query")){
    if(!arrow_is_live(cBds)){
      rlang::abort(
        "The supplied Arrow dataset has an invalid/expired C++ pointer (e.g., from a previous R session). Pass a directory path or re-open it with `arrow::open_dataset()`.",
        class = "ezbakr_arrow_disconnected"
      )
    }

    cBds <- strip_arrow_r_metadata(cBds)
    base_ds <- unwrap_arrow_dataset(cBds)
    raw_files <- tryCatch(base_ds$files, error = function(e) character(0))
    files <- if(length(raw_files) > 0L){
      normalizePath(raw_files, winslash = "/", mustWork = FALSE)
    }else{
      character(0)
    }

    fmt <- if(!is.null(format)){
      format
    }else{
      tryCatch(base_ds$format$type, error = function(e) infer_arrow_format_from_paths(files))
    }

    # Infer whether `base_ds` corresponds cleanly to a single root directory
    source_type <- "files"
    root_path <- NULL
    if(length(files) > 0L){
      has_hive <- any(grepl("/[^/=]+=[^/]+/", files))
      if(has_hive){
        cand_roots <- unique(sub("/[^/=]+=[^/]+/.*$", "", files))
        if(length(cand_roots) == 1L && dir.exists(cand_roots)){
          all_in_root <- normalizePath(
            list.files(cand_roots, recursive = TRUE, full.names = TRUE),
            winslash = "/", mustWork = FALSE
          )
          if(setequal(all_in_root, files)){
            source_type <- "directory"
            root_path <- cand_roots
          }
        }
      }else{
        cand_dirs <- unique(dirname(files))
        if(length(cand_dirs) == 1L && dir.exists(cand_dirs)){
          all_in_dir <- normalizePath(
            list.files(cand_dirs, recursive = TRUE, full.names = TRUE),
            winslash = "/", mustWork = FALSE
          )
          if(setequal(all_in_dir, files)){
            source_type <- "directory"
            root_path <- cand_dirs
          }
        }
      }
    }

    return(list(
      ds = cBds,
      source_type = source_type,
      root_path = root_path,
      explicit_paths = files,
      format = fmt,
      partitioning_cols = if(is.character(partitioning)) partitioning else NULL,
      hive_style = hive_style,
      open_dots = list()
    ))
  }

  rlang::abort(
    "`cBds` must be an Arrow `Dataset`, an `arrow_dplyr_query`, or a character path to dataset directory/files!",
    class = "cBds_invalid_type"
  )
}


# Internal helper: extract sample names from Hive partition paths or fallback query
extract_arrow_samples <- function(cBds, files){
  norm_files <- gsub("\\\\", "/", files)
  is_sample_partitioned <- length(norm_files) > 0L &&
    all(grepl("(?:^|/)sample=[^/]+(?:/|$)", norm_files))

  samps_cB <- NA_character_
  if(is_sample_partitioned && !inherits(cBds, "arrow_dplyr_query")){
    samps_cB <- unique(vapply(norm_files, function(fp){
      m <- regmatches(fp, regexec("(?:^|/)sample=([^/]+)(?:/|$)", fp))[[1]]
      if(length(m) > 1L) utils::URLdecode(m[2L]) else NA_character_
    }, character(1)))
  }

  if(anyNA(samps_cB) || length(samps_cB) == 0L){
    sample <- NULL
    samps_cB <- cBds %>%
      dplyr::select(sample) %>%
      dplyr::distinct() %>%
      dplyr::collect() %>%
      dplyr::pull(sample) %>%
      as.character()
  }

  list(
    samples = unique(as.character(samps_cB)),
    partitioned_by_sample = is_sample_partitioned
  )
}


# Internal helper: build pure-R serializable Arrow spec + mutable session cache
build_arrow_spec <- function(resolved, cB_cols, mutcounts_in_cB, basecounts_in_cB,
                             features_in_cB, samples_info){
  ds <- resolved$ds
  base_ds <- unwrap_arrow_dataset(ds)
  raw_files <- tryCatch(base_ds$files, error = function(e) character(0))
  files <- if(identical(resolved$source_type, "files") && length(resolved$explicit_paths) > 0L){
    resolved$explicit_paths
  }else if(length(raw_files) > 0L){
    normalizePath(raw_files, winslash = "/", mustWork = FALSE)
  }else{
    resolved$explicit_paths
  }

  has_unserializable_filter <- FALSE
  selected_cols <- NULL
  if(inherits(ds, "arrow_dplyr_query") && !identical(resolved$source_type, "files")){
    if(!isTRUE(ds$filtered_rows) ||
       length(ds$group_by_vars) > 0L ||
       !is.null(ds$head) ||
       !is.null(ds$tail)){
      has_unserializable_filter <- TRUE
    }
  }

  cache_env <- new.env(parent = emptyenv())
  cache_env$ds <- ds

  display_path <- if(!is.null(resolved$root_path)){
    resolved$root_path
  }else if(length(files) > 0L){
    common_d <- unique(dirname(files))
    if(length(common_d) == 1L) common_d else files[1L]
  }else{
    "<unknown>"
  }

  list(
    source_type = resolved$source_type,
    root_path = resolved$root_path,
    files = files,
    display_path = display_path,
    format = resolved$format,
    partitioning_cols = resolved$partitioning_cols,
    hive_style = isTRUE(resolved$hive_style),
    open_dots = resolved$open_dots,
    partitioned_by_sample = isTRUE(samples_info$partitioned_by_sample),
    has_unserializable_filter = has_unserializable_filter,
    cB_cols = cB_cols,
    mutcounts = mutcounts_in_cB,
    basecounts = basecounts_in_cB,
    features = features_in_cB,
    samples = samples_info$samples,
    n_files = length(files),
    cache = cache_env
  )
}


# Internal helper: open a vector of file paths while preserving Hive partition columns
open_arrow_files_source <- function(norm_paths, fmt, partitioning = NULL,
                                    hive_style = TRUE, dots = list()){
  norm_paths <- gsub("\\\\", "/", norm_paths)
  has_hive_segs <- isTRUE(hive_style) && is.null(partitioning) &&
    all(grepl("(^|/)[^/=]+=[^/]+/[^/]+$", norm_paths))

  if(has_hive_segs){
    cand_roots <- unique(sub("(?:/[^/=]+=[^/]+)+/[^/]+$", "", norm_paths))
    if(length(cand_roots) == 1L && dir.exists(cand_roots)){
      # Fast path: if file vector matches all files under cand_roots, open via root
      root_args <- list(sources = cand_roots, format = fmt, hive_style = hive_style)
      if(length(dots) > 0L) root_args <- c(root_args, dots)
      ds_root <- tryCatch(do.call(arrow::open_dataset, root_args), error = function(e) NULL)
      if(arrow_is_live(ds_root)){
        root_files <- normalizePath(ds_root$files, winslash = "/", mustWork = FALSE)
        if(setequal(root_files, norm_paths)){
          return(ds_root)
        }
        if(all(norm_paths %in% root_files) && "sample" %in% names(arrow::schema(ds_root))){
          parse_hive_sample <- function(p){
            vapply(p, function(fp){
              m <- regmatches(fp, regexec("(?:^|/)sample=([^/]+)(?:/|$)", fp))[[1]]
              if(length(m) > 1L) utils::URLdecode(m[2L]) else NA_character_
            }, character(1))
          }
          sel_samps <- unique(parse_hive_sample(norm_paths))
          root_samps <- parse_hive_sample(root_files)
          if(!anyNA(sel_samps) && setequal(root_files[root_samps %in% sel_samps], norm_paths)){
            sample <- NULL
            return(ds_root %>% dplyr::filter(sample %in% sel_samps))
          }
        }
      }
    }
  }

  open_args <- list(sources = norm_paths, format = fmt, hive_style = hive_style)
  if(!is.null(partitioning)){
    open_args$partitioning <- partitioning
  }
  if(length(dots) > 0L){
    open_args <- c(open_args, dots)
  }
  do.call(arrow::open_dataset, open_args)
}


# Internal helper: attempt to reopen an Arrow dataset from its pure-R spec
reopen_arrow_from_spec <- function(spec){
  if(is.null(spec) || isTRUE(spec$has_unserializable_filter)){
    return(NULL)
  }

  fmt <- if(!is.null(spec$format)) spec$format else "parquet"
  hive_st <- if(!is.null(spec$hive_style)) isTRUE(spec$hive_style) else TRUE
  dots <- if(is.list(spec$open_dots)) spec$open_dots else list()

  open_with_sources <- function(src){
    open_args <- list(sources = src, format = fmt, hive_style = hive_st)
    if(!is.null(spec$partitioning_cols)){
      open_args$partitioning <- spec$partitioning_cols
    }
    if(length(dots) > 0L){
      open_args <- c(open_args, dots)
    }
    tryCatch(do.call(arrow::open_dataset, open_args), error = function(e) NULL)
  }

  if(identical(spec$source_type, "directory") &&
     !is.null(spec$root_path) &&
     dir.exists(spec$root_path)){
    ds <- open_with_sources(spec$root_path)
    if(arrow_is_live(ds)) return(strip_arrow_r_metadata(ds))
  }

  if(!is.null(spec$files) && length(spec$files) > 0L && all(file.exists(spec$files))){
    if(!identical(spec$source_type, "files") && isTRUE(spec$partitioned_by_sample)){
      cand_roots <- unique(sub("/[^/=]+=[^/]+/.*$", "", gsub("\\\\", "/", spec$files)))
      if(length(cand_roots) == 1L && dir.exists(cand_roots)){
        ds_root <- open_with_sources(cand_roots)
        if(arrow_is_live(ds_root)) return(strip_arrow_r_metadata(ds_root))
      }
    }
    ds <- tryCatch(
      open_arrow_files_source(spec$files,
                              fmt = fmt,
                              partitioning = spec$partitioning_cols,
                              hive_style = hive_st,
                              dots = dots),
      error = function(e) NULL
    )
    if(arrow_is_live(ds)) return(strip_arrow_r_metadata(ds))
  }

  NULL
}


# Internal accessor: return live Arrow dataset, transparently auto-reopening after readRDS()
get_arrow_ds <- function(obj, require_live = TRUE){
  if(arrow_is_live(obj$cBds)){
    return(strip_arrow_r_metadata(obj$cBds))
  }

  spec <- obj[["metadata"]][["arrow_spec"]]
  if(!is.null(spec)){
    if(is.environment(spec$cache) && arrow_is_live(spec$cache$ds)){
      return(spec$cache$ds)
    }

    reopened <- reopen_arrow_from_spec(spec)
    if(!is.null(reopened) && arrow_is_live(reopened)){
      if(is.environment(spec$cache)){
        spec$cache$ds <- reopened
      }
      return(reopened)
    }
  }

  if(!require_live){
    return(NULL)
  }

  expected_loc <- if(!is.null(spec$display_path)){
    spec$display_path
  }else{
    "<unrecorded path>"
  }

  extra_note <- if(!is.null(spec) && isTRUE(spec$has_unserializable_filter)){
    " Because the original dataset used an in-memory `dplyr::filter()` query, its C++ filter expression cannot be auto-restored from RDS."
  }else{
    paste0(" Original dataset location not found on disk: '", expected_loc, "'.")
  }

  rlang::abort(
    paste0(
      "The Arrow dataset backing this `EZbakRArrowData` object is disconnected (likely after `readRDS()`).",
      extra_note,
      "\nReconnect it non-destructively with:\n  obj <- ReconnectArrow(obj, \"path/to/arrow_dataset\")"
    ),
    class = "ezbakr_arrow_disconnected"
  )
}


# Internal helper: verify an Arrow column holds whole numbers >= min_val
check_arrow_whole_col <- function(cBds, cBschema, col, min_val = 0){
  tname <- tolower(cBschema[[col]]$type$name)
  if(grepl("^uint", tname)){
    if(min_val <= 0) return(TRUE)
    bad <- cBds %>%
      dplyr::select(dplyr::all_of(col)) %>%
      dplyr::filter(is.na(.data[[col]]) | .data[[col]] < min_val) %>%
      utils::head(1) %>%
      dplyr::collect()
    return(nrow(bad) == 0L)
  }
  if(grepl("^int", tname)){
    if(min_val <= 0){
      bad <- cBds %>%
        dplyr::select(dplyr::all_of(col)) %>%
        dplyr::filter(is.na(.data[[col]]) | .data[[col]] < min_val) %>%
        utils::head(1) %>%
        dplyr::collect()
      return(nrow(bad) == 0L)
    }else{
      bad <- cBds %>%
        dplyr::select(dplyr::all_of(col)) %>%
        dplyr::filter(is.na(.data[[col]]) | .data[[col]] < min_val) %>%
        utils::head(1) %>%
        dplyr::collect()
      return(nrow(bad) == 0L)
    }
  }
  if(grepl("^(float|double|halffloat)", tname)){
    bad <- cBds %>%
      dplyr::select(dplyr::all_of(col)) %>%
      dplyr::filter(is.na(.data[[col]]) |
                      .data[[col]] < min_val |
                      floor(.data[[col]]) != .data[[col]]) %>%
      utils::head(1) %>%
      dplyr::collect()
    return(nrow(bad) == 0L)
  }
  FALSE
}


#' `EZbakRArrowData` object constructor for internal use
#'
#' \code{new_EZbakRArrowData} efficiently creates an object of class `EZbakRArrowData`.
#' It does not perform any rigorous checks of the legitimacy of this object.
#' @param cBds Arrow dataset (or `arrow_dplyr_query`) tracking the sample ID,
#' mutational and nucleotide content, and feature assignment of sequencing reads.
#' @param metadf Data frame tracking features of each of the samples included
#' in `cBds`.
#' @param arrow_spec Optional internal list storing serializable Arrow dataset metadata.
new_EZbakRArrowData <- function(cBds, metadf, arrow_spec = NULL){
  stopifnot((methods::is(cBds, "ArrowObject") && methods::is(cBds, "Dataset")) ||
              inherits(cBds, "arrow_dplyr_query"))
  stopifnot(is.data.frame(metadf))
  out <- list(cBds = cBds, metadf = metadf)
  if(!is.null(arrow_spec)){
    out$metadata <- list(arrow_spec = arrow_spec)
  }
  structure(out, class = "EZbakRArrowData")
}


#' `EZbakRArrowData` validator
#'
#' \code{validate_EZbakRArrowData} ensures that input for `EZbakRArrowData` object construction
#' is valid.
#' @param obj An object of class `EZbakRArrowData`
#' @param resolved Optional resolved Arrow input metadata from `resolve_arrow_input()`.
validate_EZbakRArrowData <- function(obj, resolved = NULL){

  mutcounts <- ezbakr_valid_mutcounts()

  ### Extract objects and info from objects
  vals <- unclass(obj)
  cBds <- strip_arrow_r_metadata(vals$cBds)
  obj$cBds <- cBds
  metadf <- dplyr::as_tibble(vals$metadf)

  if(is.null(resolved)){
    resolved <- resolve_arrow_input(cBds)
  }

  metadf_cols <- colnames(metadf)
  cBschema <- arrow::schema(cBds)
  cB_cols <- names(cBschema)

  # Which mutcounts are in the cB?
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]

  ### Does cB_cols contain sample and n?
  if(!("sample" %in% cB_cols & "n" %in% cB_cols)){

    rlang::abort(
      "cBds must include columns named sample (representing the sample ID)
      and n (representing the number of reads with identical values for the
      columns included in the cBds)!",
      class = "cBds_sample_n"
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
      class = "cBds_basecount"
    )

  }

  ### Make sure there is at least one feature column in cBds
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB, basecounts_expected, "sample", "n"))]
  if(length(features_in_cB) == 0L){
    rlang::abort(
      "No genomic feature columns were detected in cBds! At least one column besides sample, n, mutation counts, and nucleotide counts is required.",
      class = "cBds_no_features"
    )
  }

  ### Does metadf contain sample and correct tl(s)?
  if(!("sample" %in% metadf_cols)){
    rlang::abort(
      "metadf must include a column named sample (representing the sample ID)!",
      class = "metadf_sample"
    )
  }

  if(length(mutcounts_in_cB) > 1){

    tl_expected <- paste0("tl_", mutcounts_in_cB)
    tl_expected_p <- paste0("tpulse_", mutcounts_in_cB)
    tl_expected_c <- paste0("tchase_", mutcounts_in_cB)

    if(!all(tl_expected %in% metadf_cols)){

      if(!(all(tl_expected_p %in% metadf_cols) & all(tl_expected_c %in% metadf_cols))){

        rlang::abort(
          "Not all of the relevant label times are included in the metadf.
          For example, if your cBds has columns TC and GA, your metadf must
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

  ### Check that label times are numeric and >= 0
  is_pos_num <- function(x){
    is.numeric(x) && all(!is.na(x) & x >= 0)
  }

  if(!all(sapply(metadf[, tl_cols, drop = FALSE], is_pos_num))){

    rlang::abort("Not all columns of metadf representing label times are numbers
         greater than or equal to 0.",
                 class = "tl_pos_num")

  }

  ### Extract sample names (supports hyphenated/dotted/URL-encoded Hive names & unpartitioned datasets)
  base_ds <- unwrap_arrow_dataset(cBds)
  files <- if(identical(resolved$source_type, "files") && length(resolved$explicit_paths) > 0L){
    resolved$explicit_paths
  }else{
    tryCatch(base_ds$files, error = function(e) character(0))
  }
  samples_info <- extract_arrow_samples(cBds, files)
  samps_cB <- samples_info$samples

  if(!samples_info$partitioned_by_sample){
    message(
      "Note: Arrow dataset is not Hive-partitioned by 'sample'. ",
      "Partitioning by 'sample' (e.g., via `WriteCBArrow()` or `arrow::write_dataset(..., partitioning = 'sample')`) ",
      "is strongly recommended for faster, lower-memory per-sample queries."
    )
  }

  if(!all(as.character(metadf$sample) %in% samps_cB)){
    rlang::abort("Not all samples in metadf are present in the cB!",
                 class = "cB_metadf_samples")
  }

  ### Check that mutation counts, base counts, and n contain whole numbers
  if(!all(vapply(mutcounts_in_cB,
                 function(col) check_arrow_whole_col(cBds, cBschema, col, min_val = 0),
                 logical(1)))){

    rlang::abort("Not all columns tracking counts of mutations are integers!",
                 class = "mut_not_int")

  }

  if(!all(vapply(unique(basecounts_expected),
                 function(col) check_arrow_whole_col(cBds, cBschema, col, min_val = 0),
                 logical(1)))){

    rlang::abort("Not all columns tracking counts of nucleotides are integers!",
                 class = "basecounts_not_pos_whole")

  }

  if(!check_arrow_whole_col(cBds, cBschema, "n", min_val = 1)){

    rlang::abort("cBds column n does not contain positive integer data!",
                 class = "cB_n_int")

  }

  ### Check that mutation counts never exceed corresponding base counts
  for(idx in seq_along(mutcounts_in_cB)){
    mut_col <- mutcounts_in_cB[idx]
    base_col <- basecounts_expected[idx]
    bad_rows <- cBds %>%
      dplyr::select(dplyr::all_of(c(mut_col, base_col))) %>%
      dplyr::filter(.data[[mut_col]] > .data[[base_col]]) %>%
      utils::head(1) %>%
      dplyr::collect()
    if(nrow(bad_rows) > 0){
      rlang::abort(
        paste0("Mutation counts in column ", mut_col,
               " exceed nucleotide counts in column ", base_col, "!"),
        class = "mut_exceeds_basecount"
      )
    }
  }

  ### Build serializable Arrow spec for transparent reload after saveRDS()/readRDS()
  spec <- build_arrow_spec(
    resolved = resolved,
    cB_cols = cB_cols,
    mutcounts_in_cB = mutcounts_in_cB,
    basecounts_in_cB = unique(basecounts_expected),
    features_in_cB = features_in_cB,
    samples_info = samples_info
  )

  if(isTRUE(spec$has_unserializable_filter)){
    message(
      "Note: `cBds` includes an in-memory `arrow_dplyr_query` filter/aggregation; ",
      "if you save this object with `saveRDS()`, call `ReconnectArrow()` after `readRDS()` before running new raw-read queries."
    )
  }

  ### Convert all other metadf columns to factors
  cols_to_convert <- metadf_cols[!(metadf_cols %in% c(tl_cols, "sample"))]

  metadf <- as.data.frame(metadf)
  metadf$sample <- as.character(metadf$sample)
  metadf[cols_to_convert] <- lapply(metadf[cols_to_convert], as.factor)
  metadf <- dplyr::as_tibble(metadf)

  obj$metadf <- metadf
  if(is.null(obj$metadata)) obj$metadata <- list()
  obj$metadata$arrow_spec <- spec

  return(obj)

}


#' `EZbakRArrowData` object constructor for on-disk Arrow datasets
#'
#' \code{EZbakRArrowData} creates an object of class `EZbakRArrowData` backed by an
#' on-disk Arrow dataset and validates the provided inputs. You may pass either an
#' open `arrow::Dataset` (or `arrow_dplyr_query`) OR a character path directly to a
#' Parquet/Feather/CSV/TSV dataset directory or vector of files.
#'
#' Objects created with `EZbakRArrowData()` record serializable dataset metadata in
#' `obj$metadata$arrow_spec` so that saving with `saveRDS()` and reloading later with
#' `readRDS()` automatically reconnects to the Arrow dataset when the files remain at
#' the same path. If the dataset folder has moved, use \code{\link{ReconnectArrow}} to
#' point the object at the new location without losing any downstream analyses.
#'
#' @param cBds Either a character path (directory path or vector of file paths) OR an
#' `arrow::Dataset` / `arrow_dplyr_query` with the following fields:
#' \itemize{
#'  \item sample: Name given to particular sample from which data was collected.
#'  \item mutational counts: Integers corresponding to the number of a particular
#'  mutation seen in a sequencing read (`TC`, `GA`, `NtoA`, etc.).
#'  \item base nucleotide count: Integers corresponding to the number of instances of a
#'  particular type of nucleotide (`nT`, `nG`, `nA`, `nC`, `nU`, `nN`).
#'  \item features: One or more genomic feature columns (`XF`, `GF`, `TEC`, etc.).
#'  \item n: Positive integer count of reads with identical values for all other columns.
#' }
#' @param metadf Data frame detailing various aspects of each of the samples included
#' in `cBds` (`sample`, `tl` / `tl_<muttype>` / `tpulse` + `tchase`, and experimental
#' covariates such as `treatment` or `batch`).
#' @param format Optional file format string passed to `arrow::open_dataset()` when
#' `cBds` is a character path (e.g., `"parquet"`, `"feather"`, `"csv"`, `"tsv"`).
#' Inferred automatically from file extensions when `NULL` (default).
#' @param partitioning Optional partitioning specification passed to
#' `arrow::open_dataset()` when `cBds` is a character path. Leave `NULL` (default) to
#' automatically detect Hive-style `sample=...` directory partitioning.
#' @param hive_style Logical passed to `arrow::open_dataset()` when `cBds` is a
#' character path (`TRUE` by default).
#' @param ... Additional arguments forwarded to `arrow::open_dataset()` when `cBds`
#' is a character path.
#' @return An `EZbakRArrowData` object containing `cBds`, `metadf`, and serializable
#' Arrow connection metadata under `metadata$arrow_spec`.
#' @import data.table
#' @examples
#' library(arrow)
#' simdata <- EZSimulate(30)
#'
#' dataset_dir <- file.path(tempdir(), "ezbakr_arrow_example")
#' WriteCBArrow(simdata$cB, path = dataset_dir, overwrite = TRUE)
#'
#' # Pass the directory path directly (or pass open_dataset(dataset_dir))
#' ezbdo <- EZbakRArrowData(dataset_dir, simdata$metadf)
#'
#' @export
EZbakRArrowData <- function(cBds,
                            metadf,
                            format = NULL,
                            partitioning = NULL,
                            hive_style = TRUE,
                            ...){

  resolved <- resolve_arrow_input(
    cBds = cBds,
    format = format,
    partitioning = partitioning,
    hive_style = hive_style,
    ...
  )

  metadf_copy <- data.table::as.data.table(metadf)

  validate_EZbakRArrowData(
    new_EZbakRArrowData(resolved$ds, metadf_copy),
    resolved = resolved
  )

}


#' Reconnect an `EZbakRArrowData` object to its on-disk Arrow dataset
#'
#' After saving an `EZbakRArrowData` object with `saveRDS()` and loading it in a new R
#' session with `readRDS()`, EZbakR automatically reconnects to the original Arrow
#' dataset path if it still exists on disk. If you moved the dataset directory (for
#' example, after transferring files from an HPC cluster to a local machine) or want
#' to attach a new `arrow::Dataset` / `arrow_dplyr_query`, call `ReconnectArrow()` to
#' re-bind the Arrow backend while preserving all existing downstream tables
#' (`fractions`, `kinetics`, `averages`, `comparisons`, `dynamics`, `readcounts`, etc.).
#'
#' @param obj An `EZbakRArrowData` object.
#' @param cBds Optional new dataset path (directory or vector of files) or live
#' `arrow::Dataset` / `arrow_dplyr_query`. If `NULL` (the default), attempts to
#' reopen the path recorded in `obj$metadata$arrow_spec`.
#' @param format Optional file format passed to `arrow::open_dataset()` when `cBds`
#' is a character path.
#' @param partitioning Optional partitioning specification passed to
#' `arrow::open_dataset()` when `cBds` is a character path.
#' @param hive_style Logical passed to `arrow::open_dataset()` (`TRUE` by default).
#' @param validate Logical; if `TRUE` (default), verify that required columns and
#' `metadf$sample` entries exist in the reconnected dataset.
#' @param ... Additional arguments forwarded to `arrow::open_dataset()` when `cBds`
#' is a character path.
#' @return The input `obj` with a live `cBds` handle and updated `metadata$arrow_spec`,
#' with all existing analysis results intact.
#' @examples
#' simdata <- EZSimulate(20)
#' ds_dir <- file.path(tempdir(), "ezbakr_reconnect_demo")
#' WriteCBArrow(simdata$cB, ds_dir, overwrite = TRUE)
#' ezbdo <- EZbakRArrowData(ds_dir, simdata$metadf)
#'
#' rds_path <- tempfile(fileext = ".rds")
#' saveRDS(ezbdo, rds_path)
#' ezbdo_loaded <- readRDS(rds_path)
#' ezbdo_loaded <- ReconnectArrow(ezbdo_loaded, ds_dir)
#'
#' @export
ReconnectArrow <- function(obj,
                           cBds = NULL,
                           format = NULL,
                           partitioning = NULL,
                           hive_style = TRUE,
                           validate = TRUE,
                           ...){

  if(!methods::is(obj, "EZbakRArrowData")){
    rlang::abort(
      "`obj` must inherit from `EZbakRArrowData`!",
      class = "reconnect_not_arrow"
    )
  }

  old_spec <- obj[["metadata"]][["arrow_spec"]]

  if(is.null(cBds)){
    reopened <- reopen_arrow_from_spec(old_spec)
    if(is.null(reopened) || !arrow_is_live(reopened)){
      loc <- if(!is.null(old_spec$display_path)) old_spec$display_path else "<unknown>"
      rlang::abort(
        paste0("Could not automatically reopen Arrow dataset at '", loc,
               "'. Please supply the new directory/file path to `ReconnectArrow(obj, 'new/path')`."),
        class = "ezbakr_arrow_disconnected"
      )
    }
    resolved <- resolve_arrow_input(
      reopened,
      format = if(!is.null(format)) format else old_spec$format,
      partitioning = if(!is.null(partitioning)) partitioning else old_spec$partitioning_cols,
      hive_style = hive_style
    )
  }else{
    resolved <- resolve_arrow_input(
      cBds = cBds,
      format = format,
      partitioning = partitioning,
      hive_style = hive_style,
      ...
    )
  }

  ds <- resolved$ds
  cBschema <- arrow::schema(ds)
  cB_cols <- names(cBschema)
  mutcounts <- ezbakr_valid_mutcounts()
  mutcounts_in_cB <- cB_cols[cB_cols %in% mutcounts]
  basecounts_expected <- paste0("n", substr(mutcounts_in_cB, start = 1, stop = 1))
  features_in_cB <- cB_cols[!(cB_cols %in% c(mutcounts_in_cB, basecounts_expected, "sample", "n"))]

  if(!all(c("sample", "n") %in% cB_cols)){
    rlang::abort(
      "Reconnected Arrow dataset must include `sample` and `n` columns!",
      class = "cBds_sample_n"
    )
  }
  if(!all(basecounts_expected %in% cB_cols)){
    rlang::abort(
      "Reconnected Arrow dataset is missing required nucleotide base-count columns!",
      class = "cBds_basecount"
    )
  }
  if(length(features_in_cB) == 0L){
    rlang::abort(
      "Reconnected Arrow dataset has no genomic feature columns!",
      class = "cBds_no_features"
    )
  }

  base_ds <- unwrap_arrow_dataset(ds)
  files <- if(identical(resolved$source_type, "files") && length(resolved$explicit_paths) > 0L){
    resolved$explicit_paths
  }else{
    tryCatch(base_ds$files, error = function(e) character(0))
  }
  samples_info <- extract_arrow_samples(ds, files)

  if(validate){
    if(!all(as.character(obj$metadf$sample) %in% samples_info$samples)){
      rlang::abort(
        "Not all samples in `obj$metadf` are present in the reconnected Arrow dataset!",
        class = "cB_metadf_samples"
      )
    }
  }

  new_spec <- build_arrow_spec(
    resolved = resolved,
    cB_cols = cB_cols,
    mutcounts_in_cB = mutcounts_in_cB,
    basecounts_in_cB = unique(basecounts_expected),
    features_in_cB = features_in_cB,
    samples_info = samples_info
  )

  obj$cBds <- ds
  if(is.null(obj$metadata)) obj$metadata <- list()
  obj$metadata$arrow_spec <- new_spec

  return(obj)

}


#' Write `cB` data to a `sample`-partitioned Parquet Arrow dataset
#'
#' Convenience helper that writes an in-memory `cB` table, an `EZbakRData` object,
#' an `arrow::Dataset`, or one or more `.csv` / `.csv.gz` / `.tsv` / `.tsv.gz` /
#' `.parquet` `cB` files (e.g., from `fastq2EZbakR`) to an on-disk Arrow dataset
#' partitioned by `sample`, automatically casting mutation counts, nucleotide base
#' counts, and `n` to compact integer types (`int32`, or `int64` when values exceed
#' 32-bit integer range).
#'
#' @param cB An in-memory `data.frame` / `tibble` / `data.table`, an `EZbakRData`
#' object, an `arrow::Dataset` / `arrow_dplyr_query`, or a character vector of
#' existing `cB` file paths (`.csv`, `.csv.gz`, `.tsv`, `.tsv.gz`, `.parquet`, etc.).
#' @param path Target directory path where the partitioned Arrow dataset will be written.
#' @param format Output Arrow file format (`"parquet"` by default).
#' @param partitioning Character vector of columns to Hive-partition by (`"sample"`
#' by default).
#' @param overwrite Logical; if `FALSE` (default) and `path` already exists and is
#' non-empty, an error is thrown to prevent mixing shards from separate runs. If
#' `TRUE`, existing contents of `path` are removed before writing.
#' @param ... Additional arguments passed to `arrow::write_dataset()`.
#' @return Invisibly returns the normalized `path` to the written dataset directory,
#' suitable for piping directly into \code{\link{EZbakRArrowData}}.
#' @examples
#' simdata <- EZSimulate(20)
#' out_dir <- file.path(tempdir(), "ezbakr_write_cb_demo")
#' WriteCBArrow(simdata$cB, out_dir, overwrite = TRUE)
#' ezbdo <- EZbakRArrowData(out_dir, simdata$metadf)
#'
#' @export
WriteCBArrow <- function(cB,
                         path,
                         format = "parquet",
                         partitioning = "sample",
                         overwrite = FALSE,
                         ...){

  sample <- NULL

  if(!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)){
    rlang::abort("`path` must be a single non-empty directory path string!")
  }

  if(dir.exists(path)){
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if(length(existing) > 0L){
      if(!overwrite){
        rlang::abort(
          paste0("Target directory '", path, "' already exists and is non-empty. ",
                 "Set `overwrite = TRUE` to replace it or choose a new directory."),
          class = "arrow_dir_exists"
        )
      }
      unlink(path, recursive = TRUE, force = TRUE)
    }
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  if(methods::is(cB, "EZbakRData")){
    cB <- cB$cB
  }else if(methods::is(cB, "EZbakRArrowData")){
    cB <- get_arrow_ds(cB)
  }

  mutcounts <- ezbakr_valid_mutcounts()

  if(is.data.frame(cB)){
    dt <- data.table::as.data.table(data.table::copy(cB))
    cols <- colnames(dt)
    if(!all(c("sample", "n") %in% cols)){
      rlang::abort("`cB` must contain columns `sample` and `n`!", class = "cBds_sample_n")
    }
    mut_in <- cols[cols %in% mutcounts]
    base_in <- intersect(paste0("n", substr(mut_in, 1, 1)), cols)
    int_cols <- unique(c(mut_in, base_in, "n"))

    dt[, sample := as.character(sample)]
    for(ic in int_cols){
      v <- dt[[ic]]
      if(!is.numeric(v) || anyNA(v) || any(v < 0) || any(floor(v) != v)){
        rlang::abort(
          paste0("Column `", ic, "` must contain non-negative whole numbers!"),
          class = "cB_non_integer"
        )
      }
      if(ic == "n" && any(v < 1)){
        rlang::abort("Column `n` must contain positive whole numbers (>= 1)!", class = "cB_n_int")
      }
      if(max(v, 0) <= .Machine$integer.max){
        data.table::set(dt, j = ic, value = as.integer(v))
      }
    }

    arrow::write_dataset(
      dataset = as.data.frame(dt),
      path = path,
      format = format,
      partitioning = partitioning,
      ...
    )
  }else if(is.character(cB) || methods::is(cB, "Dataset") || inherits(cB, "arrow_dplyr_query")){
    resolved <- resolve_arrow_input(cB)
    ds_in <- resolved$ds
    sc <- arrow::schema(ds_in)
    cols <- names(sc)
    if(!all(c("sample", "n") %in% cols)){
      rlang::abort("Input dataset/files must contain columns `sample` and `n`!", class = "cBds_sample_n")
    }
    mut_in <- cols[cols %in% mutcounts]
    base_in <- intersect(paste0("n", substr(mut_in, 1, 1)), cols)
    int_cols <- unique(c(mut_in, base_in, "n"))

    q_cast <- dplyr::mutate(
      ds_in,
      sample = as.character(sample),
      dplyr::across(dplyr::all_of(int_cols), as.integer)
    )

    arrow::write_dataset(
      dataset = q_cast,
      path = path,
      format = format,
      partitioning = partitioning,
      ...
    )
  }else{
    rlang::abort(
      "`cB` must be a data.frame/data.table, `EZbakRData`, `arrow::Dataset`, or character vector of file paths!"
    )
  }

  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))

}


#' Materialize an `EZbakRArrowData` object (or filtered slice) into an in-memory `EZbakRData` object
#'
#' Pulls read-level `cB` data from an on-disk `EZbakRArrowData` object into a standard
#' in-memory \code{\link{EZbakRData}} object, pushing optional `samples` and `features`
#' filters down into Arrow before collecting so you can inspect raw mutation counts or
#' likelihood curves for specific genes/samples without loading the entire dataset into RAM.
#'
#' @param obj An `EZbakRArrowData` object.
#' @param samples Optional character vector of sample names to collect. Defaults to
#' all samples in `obj$metadf$sample`.
#' @param features Optional character vector of feature IDs (e.g., gene symbols or
#' transcript IDs) to filter to before collecting into R memory.
#' @param feature_col Optional name of the feature column in `cBds` to filter against
#' when `features` is supplied. Required only if `cBds` contains multiple feature columns.
#' @param copy_analyses Logical; if `TRUE` (default), copy existing `fractions`,
#' `kinetics`, `mutation_rates`, `readcounts`, and when unfiltered (`samples = NULL`
#' and `features = NULL`) `averages`, `comparisons`, and `dynamics` slots onto the
#' returned `EZbakRData` object.
#' @return An in-memory `EZbakRData` object.
#' @examples
#' simdata <- EZSimulate(25)
#' ds_dir <- file.path(tempdir(), "ezbakr_collect_demo")
#' WriteCBArrow(simdata$cB, ds_dir, overwrite = TRUE)
#' ez_arrow <- EZbakRArrowData(ds_dir, simdata$metadf)
#'
#' # Materialize two features for interactive inspection
#' ez_sub <- CollectArrow(ez_arrow, features = c("feature1", "feature2"))
#'
#' @export
CollectArrow <- function(obj,
                         samples = NULL,
                         features = NULL,
                         feature_col = NULL,
                         copy_analyses = TRUE){

  sample <- NULL

  if(!methods::is(obj, "EZbakRArrowData")){
    rlang::abort("`obj` must be an `EZbakRArrowData` object!", class = "collect_not_arrow")
  }

  cBds <- get_arrow_ds(obj, require_live = TRUE)
  cBschema <- arrow::schema(cBds)
  cB_cols <- names(cBschema)
  mutcounts <- ezbakr_valid_mutcounts()
  mut_in <- cB_cols[cB_cols %in% mutcounts]
  base_in <- paste0("n", substr(mut_in, 1, 1))
  feat_cols <- cB_cols[!(cB_cols %in% c(mut_in, base_in, "sample", "n"))]

  metadf_sub <- dplyr::as_tibble(obj$metadf)
  if(!is.null(samples)){
    samples <- as.character(samples)
    if(!all(samples %in% as.character(metadf_sub$sample))){
      rlang::abort("Some requested `samples` are not present in `obj$metadf$sample`!")
    }
    metadf_sub <- metadf_sub[as.character(metadf_sub$sample) %in% samples, , drop = FALSE]
  }
  keep_samples <- as.character(metadf_sub$sample)

  q <- cBds %>% dplyr::filter(sample %in% keep_samples)

  if(!is.null(features)){
    if(is.null(feature_col)){
      if(length(feat_cols) == 1L){
        feature_col <- feat_cols[1L]
      }else{
        rlang::abort(
          paste0("Multiple feature columns exist in `cBds` (",
                 paste(feat_cols, collapse = ", "),
                 "); please specify `feature_col` when filtering by `features`.")
        )
      }
    }else if(!(feature_col %in% feat_cols)){
      rlang::abort(paste0("`feature_col` ('", feature_col, "') is not a feature column in `cBds`!"))
    }
    features_chr <- as.character(features)
    q <- q %>% dplyr::filter(.data[[feature_col]] %in% features_chr)
  }

  cB_mem <- dplyr::collect(q)
  if(nrow(cB_mem) == 0L){
    rlang::abort("No rows in the Arrow dataset matched the requested `samples` / `features` filter!")
  }

  # Ensure integer-compatible numeric columns and plain tibble
  cB_mem <- dplyr::as_tibble(cB_mem)
  cB_mem$sample <- as.character(cB_mem$sample)

  # Ensure metadf only keeps samples actually present when feature filtering drops a sample completely
  present_samps <- unique(cB_mem$sample)
  metadf_sub <- metadf_sub[as.character(metadf_sub$sample) %in% present_samps, , drop = FALSE]

  ez_mem <- EZbakRData(cB_mem, metadf_sub)

  if(copy_analyses){
    unfiltered <- is.null(samples) && is.null(features)

    filter_sample_feat_tbl <- function(tbl){
      out <- tbl
      if("sample" %in% colnames(out)){
        out <- out[as.character(out$sample) %in% present_samps, , drop = FALSE]
      }
      if(!is.null(features) && !is.null(feature_col) && (feature_col %in% colnames(out))){
        out <- out[as.character(out[[feature_col]]) %in% as.character(features), , drop = FALSE]
      }
      out
    }

    for(slot_nm in c("mutation_rates", "fractions", "kinetics", "readcounts", "dropout_params")){
      if(!is.null(obj[[slot_nm]])){
        ez_mem[[slot_nm]] <- lapply(obj[[slot_nm]], function(el){
          if(is.data.frame(el)) filter_sample_feat_tbl(el) else el
        })
        if(!is.null(obj$metadata[[slot_nm]])){
          ez_mem$metadata[[slot_nm]] <- obj$metadata[[slot_nm]]
        }
      }
    }

    if(unfiltered){
      for(slot_nm in c("averages", "comparisons", "dynamics")){
        if(!is.null(obj[[slot_nm]])){
          ez_mem[[slot_nm]] <- obj[[slot_nm]]
          if(!is.null(obj$metadata[[slot_nm]])){
            ez_mem$metadata[[slot_nm]] <- obj$metadata[[slot_nm]]
          }
        }
      }
    }else if(!is.null(features) && is.null(samples)){
      for(slot_nm in c("averages", "comparisons", "dynamics")){
        if(!is.null(obj[[slot_nm]])){
          ez_mem[[slot_nm]] <- lapply(obj[[slot_nm]], function(el){
            if(is.data.frame(el) && !is.null(feature_col) && (feature_col %in% colnames(el))){
              el[as.character(el[[feature_col]]) %in% as.character(features), , drop = FALSE]
            }else{
              el
            }
          })
          if(!is.null(obj$metadata[[slot_nm]])){
            ez_mem$metadata[[slot_nm]] <- obj$metadata[[slot_nm]]
          }
        }
      }
    }

    extra_classes <- setdiff(class(obj), c("EZbakRArrowData", "EZbakRData"))
    if(!unfiltered && !is.null(samples)){
      extra_classes <- setdiff(extra_classes, c("EZbakRFit", "EZbakRCompare", "EZbakRDynamics"))
    }
    class(ez_mem) <- unique(c(extra_classes, "EZbakRData"))
  }

  return(ez_mem)

}
