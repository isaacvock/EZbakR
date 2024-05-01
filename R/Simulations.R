SimulateOneRep <- function(nfeatures, read_vect = NULL, label_time = 2,
                           sample_name = "sampleA",
                           feature_prefix = "Gene",
                           fn_vect = NULL, kdeg_vect = NULL, ksyn_vect = NULL,
                           pnew = 0.05, pold = 0.002,
                           logkdeg_mean = -1.9, logkdeg_sd = 0.7,
                           logksyn_mean = 2.3, logksyn_sd = 0.7,
                           seqdepth = 10000000){
  
  ### Check validity of input
  
  args <- c(as.list(environment()))
  
  check_SimulateOneRep_input(args)
  
  
  
  
  ### Simulate kinetic parameters as necesary
  
  # kdeg and fraction new
  
  if(is.null(fn_vect)){
    
    if(is.null(kdeg_vect)){
      
      kdeg_vect <- rlnorm(nfeatures,
                          logkdeg_mean,
                          logkdeg_sd)
      
      
    }
    
    fn_vect <- 1 - exp(-kdeg_vect*label_time)
    
    
  }
  
  # read counts
  
  if(is.null(read_vect)){
    
    if(is.null(ksyn_vect)){
      
      ksyn_vect <- rlnorm(nfeatures,
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
  
  read_status <- rbinom(n = totreads,
                        size = 1,
                        prob = rep(fn_vect, times = read_vect))
  
  nT_count <- rbinom(n = totreads,
                     size = readlength,
                     prob = Ucont)
  
  TC_count <- rbinom(n = totreads,
                     size = nT_count,
                     prob = read_status*pnew + (1 - read_status)*pold)
  
  cB <- data.table(
    sample = sample_name,
    feature = rep(paste0(feature_prefix, 1:nfeatures), 
                  times = read_vect),
    TC = TC_count,
    nT = nT_count
  )[,.(n = .N), by = .(sample, feature, TC, nT)]
  
  
  ### Save ground truth
  
  truth <- data.table(sample = sample_name,
                      feature = paste0(feature_prefix, 1:nfeatures),
                      true_fraction_highTC = fn_vect,
                      true_kdeg = kdeg_vect,
                      true_ksyn = ksyn_vect)
  
  
  return(list(cB = cB,
              ground_truth = truth))
  
  
  
}




SimulateMultiCondition <- function(nfeatures, metadf, 
                                   param_details, mean_formula, 
                                   dispslope = 5, dispint = 0.01,
                                   logkdegsdtrend_slope = -0.3,
                                   logkdegsdtrend_intercept = -1.5,
                                   logksynsdtrend_slope = -0.3,
                                   logksynsdtrend_intercept = -1.5,
                                   seqdepth = 10000000, label_time = 2,
                                   pnew = 0.05, pold = 0.001,
                                   feature_prefix = "Gene"){
  
  
  
  ### Fill metadf with parameters that are only specified as single value
  
  mcols <- colnames(metadf)
  
  if(!('seqdepth' %in% mcols)){
    
    metadf$seqdepth = seqdepth
    
  }
  
  if(!('label_time' %in% mcols)){
    
    metadf$label_time = label_time
    
  }
  
  if(!('pnew' %in% mcols)){
    
    metadf$pnew = pnew
    
  }
  
  
  if(!('pold' %in% mcols)){
    
    metadf$pold = pold
    
  }
  
  
  
  
  
  
  ### Need to simulate linear model parameter values for all parameters specified
  
  
  mean_design <- model.matrix(mean_formula, metadf)
  
  mean_design_cols <- colnames(mean_design)
  
  
  # Reference log(kdegs)
  pdref <- param_details %>%
    filter(reference)
  
  logkdeg_ref <- rnorm(nfeatures,
                       mean = pdref$logkdeg_mean,
                       sd = pdref$logkdeg_sd)
  
  
  # Reference log(ksyns)
  logksyn_ref <- rnorm(nfeatures,
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
      filter(param == mean_design_cols[p])
    
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
      
      logkdeg_params[[p]] <- is_kdeg_param_nonzero*rnorm(nfeatures,
                                                         pd$logkdeg_mean,
                                                         pd$logkdeg_sd) + 
        logkdeg_ref
      
      logksyn_params[[p]] <- is_ksyn_param_nonzero*rnorm(nfeatures,
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
                               n = nfeatures) {
    abundance <- vector("list", n) 
    logkdeg <- abundance
    logksyn <- abundance 
    
    for (i in 1:n) {
      logksyn_i <- X %*% extract_ith(logksyn_params, i)  
      logkdeg_i <- X %*% extract_ith(logkdeg_params, i)  
      
      logkdeg[[i]] <- logkdeg_i 
      logksyn[[i]] <- logksyn_i
      abundance[[i]] <- exp(logksyn_i)/exp(logkdeg_i) 
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
      means[i] <- mean(abundances_i)
      sds[i] <- sd(abundances_i)
      
    }
    
    readcounts <- lapply(abundances,
                         function(vector) (vector/sums)*seqdepths)
    zscores <- lapply(abundances,
                      function(vector) (vector - means)/sds)
    return(list(readcounts = readcounts,
                read_zscores = zscores))
    
    
  }
  
  
  ### Determine sample-specific averages to simulate from
  
  kinetics <- compute_kinetics(mean_design, logkdeg_params,
                               logksyn_params)
  
  
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
                         logkdegsdtrend_slope)
    
    logksyn_sds <- exp(kinetics_and_reads$read_zscores[[i]]*logksynsdtrend_slope +
                         logksynsdtrend_slope)
    
    
    logkdegs[[i]] <- rnorm(n = nsamp,
                           mean = feature_logkdegs,
                           sd = logkdeg_sds)
    
    logksyns[[i]] <- rnorm(n = nsamp,
                           mean = feature_logksyns,
                           sd = logksyn_sds)
    
    reads[[i]] <- rnbinom(n = nsamp,
                          mu = feature_readavgs,
                          size = 1/(dispslope/feature_readavgs + dispint))
    
  }
  
  
  ### Simulate data for each replicate
  simdata <- vector(mode = "list", length = nrow(metadf))
  for(s in 1:nrow(metadf)){
    
    simdata[[s]] <- SimulateOneRep(nfeatures = nfeatures,
                                   read_vect = extract_ith(reads, s),
                                   label_time = as.numeric(metadf[s,"label_time"]),
                                   sample_name = as.character(metadf[s, "sample"]),
                                   kdeg_vect = exp(extract_ith(logkdegs, s)),
                                   ksyn_vect = exp(extract_ith(logksyns, s)),
                                   pnew = metadf[s, "pnew"],
                                   pold = metadf[s, "pold"])
    
  }
  
  # Combine replicate simulation data objects into one 
  names_to_bind <- names(simdata[[1]])
  
  final_simdata <- lapply(names_to_bind, function(name) {
    bind_rows(lapply(simdata, function(inner_list) inner_list[[name]]))
  })
  
  names(final_simdata) <- names_to_bind
  
  
  ### Save linear model parameters in convenient format
  
  names(logkdeg_params) <- paste0("true_logkdeg_", mean_design_cols)
  names(logksyn_params) <- paste0("true_logksyn_", mean_design_cols)
  
  feature_prefix <- "Gene"
  
  kinetic_parameters <- bind_cols(list(logkdeg_params, logksyn_params))
  kinetic_parameters[['feature']] <- paste0(feature_prefix, 1:nfeatures)
  
  
  
  ### Gather output
  
  output <- list(cB = final_simdata$cB,
                 PerRepTruth = final_simdata$ground_truth,
                 AvgTruth = kinetic_parameters,
                 metadf = metadf,
                 param_details = param_details)
  
  return(output)
  
  
}






###### PARAMETER CHECKS

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
    
    stop("logkdeg_mean must be numeric")
    
  }
  
  if(lkd_sd <= 0){
    
    stop("logkdeg_sd must be >= 0!")
    
  }
  
  
  ### logksyn_mean
  lks <- args$logksyn_mean
  
  if(!is.numeric(lks)){
    
    stop("logksyn_mean must be numeric")
    
  }
  
  
  ### logkdeg_sd
  lks_sd <- args$logksyn_sd
  
  if(!is.numeric(lks_sd)){
    
    stop("logksyn_mean must be numeric")
    
  }
  
  if(lks_sd <= 0){
    
    stop("logksyn_sd must be >= 0!")
    
  }
  
  
  ### seqdepth
  sdep <- args$seqdepth
  
  if(!is.numeric(sdep)){
    
    stop("seqdepth must be numeric")
    
  }
  
  if(sdep <= 0){
    
    stop("seqdepth must be > 0!")
    
  }
  
  if(round(sdep) != sdep){
    
    stop("seqdepth must be an integer!")
    
  }
  
  
}