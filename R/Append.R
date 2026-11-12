EZappend <- function(obj, new_table,
                     features = NULL,
                     name = NULL,
                     suffix = "_obj2",
                     table_type = c("fractions", "kinetics"),
                     bind_cBs = FALSE){

  table_type <- match.arg(table_type)

  if(table_type == "fractions"){

    obj <- EZconcat(obj, EZbakRFractions(new_table, obj[['metadf']],
                                         name = name),
                    suffix = suffix,
                    bind_cBs = bind_cBs)


  }else{

    if(is.null(features)){

      stop("features must be specified if table_type is 'kinetics'")

    }

    obj <- EZconcat(obj, EZbakRKinetics(new_table, obj[['metadf']],
                                        features = features,
                                        name = name),
                    suffix = suffix,
                    bind_cBs = bind_cBs)

  }

  return(obj)

}


EZconcat <- function(obj1, obj2,
                     suffix = "_obj2",
                     bind_cBs = FALSE){

  ### What objects are in each?

  names_1 <- names(obj1)
  names_2 <- names(obj2)

  obj_to_concat <- intersect(names_1, names_2)

  obj_unique_to_1 <- setdiff(names_1[!(names_1 %in% obj_to_concat)], c("metadf", "metadata"))
  obj_unique_to_2 <- setdiff(names_2[!(names_2 %in% obj_to_concat)], c("metadf", "metadata"))

  obj_to_concat <- obj_to_concat[!(obj_to_concat %in% c("metadf", "metadata"))]


  ### Make sure metadfs are identical

  m1 <- obj1[['metadf']]
  m2 <- obj2[['metadf']]

  if(!setequal(colnames(m1), colnames(m2))){
    stop("metadfs in obj1 and obj2 are not identical!")
  }

  compare_metadf_1 <- dplyr::anti_join(m1, m2, by = colnames(m1))
  compare_metadf_2 <- dplyr::anti_join(m2, m1, by = colnames(m1))

  if(nrow(compare_metadf_1) > 0 || nrow(compare_metadf_2) > 0){
    stop("metadfs in obj1 and obj2 are not identical!")
  }


  ### Loop through each list to concatenate

  final_obj <- list(metadf = m1)

  for(name in obj_to_concat){

    l1 <- obj1[[name]]
    l2 <- obj2[[name]]

    if(name == "cB"){

      if(bind_cBs){

        final_obj[[name]] <- dplyr::bind_rows(l1, l2)

      }else{

        final_obj[[name]] <- l1

      }



    }else if(name == "cBds"){

      final_obj[[name]] <- l1
      if(!is.null(obj1[['metadata']][['arrow_spec']])){
        final_obj[['metadata']][['arrow_spec']] <- obj1[['metadata']][['arrow_spec']]
      }else if(!is.null(obj2[['metadata']][['arrow_spec']])){
        final_obj[['metadata']][['arrow_spec']] <- obj2[['metadata']][['arrow_spec']]
      }

    }else{



      ml1 <- obj1[['metadata']][[name]]
      ml2 <- obj2[['metadata']][[name]]


      names_l1 <- names(l1)
      names_l2 <- names(l2)

      ### Will preseve all names of obj1 and will try to preserve all names
      ### of obj2. If names in obj2 match that of a name in obj1, the string
      ### suffix will be appended

      newnames <- ifelse(names_l2 %in% names_l1,
                         paste0(names_l2, suffix),
                         names_l2)

      names(l2) <- newnames

      if(!is.null(ml2)){
        names(ml2) <- newnames
      }


      final_obj[[name]] <- c(l1, l2)
      if(!is.null(ml1) || !is.null(ml2)){
        final_obj[['metadata']][[name]] <- c(ml1, ml2)
      }

    }



  }

  final_obj[['metadf']] <- m1

  for(name in obj_unique_to_1){

    final_obj[[name]] <- obj1[[name]]
    if(name == "cBds" && !is.null(obj1[['metadata']][['arrow_spec']])){
      final_obj[['metadata']][['arrow_spec']] <- obj1[['metadata']][['arrow_spec']]
    }
    if(!is.null(obj1[['metadata']][[name]])){
      final_obj[['metadata']][[name]] <- obj1[['metadata']][[name]]
    }

  }

  for(name in obj_unique_to_2){

    final_obj[[name]] <- obj2[[name]]
    if(name == "cBds" && !is.null(obj2[['metadata']][['arrow_spec']])){
      final_obj[['metadata']][['arrow_spec']] <- obj2[['metadata']][['arrow_spec']]
    }
    if(!is.null(obj2[['metadata']][[name]])){
      final_obj[['metadata']][[name]] <- obj2[['metadata']][[name]]
    }

  }

  combined_classes <- unique(c(class(obj1), class(obj2)))
  canonical_order <- c("EZbakRGoF", "EZbakRDynamics", "EZbakRCompare",
                       "EZbakRFit", "EZbakRKinetics", "EZbakRFractions",
                       "EZbakRArrowData", "EZbakRData")
  ordered_classes <- intersect(canonical_order, combined_classes)
  if(length(ordered_classes) > 0){
    class(final_obj) <- ordered_classes
  }else{
    class(final_obj) <- class(obj1)
  }

  return(final_obj)

}
