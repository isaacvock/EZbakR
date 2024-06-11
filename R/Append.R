EZappend <- function(obj, new_table,
                     features = NULL,
                     name = NULL,
                     suffix = "_obj2",
                     table_type = c("fractions", "kinetics")){

  table_type <- match.arg(table_type)

  if(table_type == "fractions"){

    obj <- EZconcat(obj, EZbakRFractions(new_table, obj[['metadf']],
                                         name = name),
                    suffix = suffix)


  }else{

    if(is.null(features)){

      stop("features must be specified if table_type is 'kinetics'")

    }

    obj <- EZconcat(obj, EZbakRKinetics(new_table, obj[['metadf']],
                                        features = features,
                                        name = name),
                    suffix = suffix)

  }

}


EZconcat <- function(obj1, obj2,
                     suffix = "_obj2",
                     bind_cBs){

  ### What objects are in each?

  names_1 <- names(obj1)
  names_2 <- names(obj2)

  obj_to_concat <- intersect(names_1, names_2)

  obj_unique_to_1 <- names_1[!(names_1 %in% obj_to_concat)]
  obj_unique_to_2 <- names_2[!(names_2 %in% obj_to_concat)]

  obj_to_concat <- obj_to_concat[!(obj_to_concat %in% c("metadf", "metadata"))]


  ### Make sure metadfs are identical

  m1 <- obj1[['metadf']]
  m2 <- obj2[['metadf']]

  compare_metadf <- dplyr::anti_join(m1, m2)

  if(nrow(compare_metadf) > 0){
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



    }else{



      ml1 <- obj1[['metadata']][[name]]
      ml2 <- obj2[['metadata']][[name]]


      names_l1 <- names(l1)
      names_l2 <- names(l2)

      ### Will preseve all names of obj1 and will try to preserve all names
      ### of obj2. If names in obj2 match that of a name in obj2, the string
      ### suffix will be appended

      newnames <- ifelse(names_l2 %in% names_l1,
                         paste0(names_l2, suffix),
                         names_l2)

      names(l2) <- newnames

      names(ml2) <- newnames


      final_obj[[name]] <- c(l1, l2)
      final_obj[['metadata']][[name]] <- c(ml1, ml2)

    }



  }

  final_obj[['metadf']] <- m1

  for(name in obj_unique_to_1){

    final_obj[[name]] <- obj1[[name]]

  }

  for(name in obj_unique_to_2){

    final_obj[[name]] <- obj1[[name]]

  }



}
