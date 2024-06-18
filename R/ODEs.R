EZfitode <- function(obj,
                     graph,
                     grouping_features,
                     sub_features,
                     parameter_names = paste0("k", 1:(nrow(graph)-1))){

  ##### ORDER OF OPERATIONS
  # 1) Infer homogeneous ODE system matrix representation (A)
  # 2) Fit model, lol.

  ##### Challenges to overcome
  # 1) Need to generalize steady-state estimation which involves figuring
  # out what the vector of zeroth-order parameters looks like
  # 2) Going to use averages as input, so need to figure out how to best
  # model average of normalized read counts.
  # 3) Still struggling with a complete, efficient generalization of inference
  # of the matrix A. NOT ANY MORE: Create a diagonal matrix from the row
  # sums of the reduced adjacency matrix and add this to the transpose of the reduced
  # adjacency matrix.
  # 4)




}
