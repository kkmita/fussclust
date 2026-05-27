#' Creates DHE (stands for "distances horizontally exploded") and DVE
#' (stands for "distances vertically exploded") matrices.
#'
#' @param A Matrix of size N x c.
#' 
#' @param vertical Boolean switch.
#' If `TRUE`, create DVE (vertical explosion).
#' If `FALSE`, create DHE (horizontal explosion).
#'
#' @return Matrix of size Nc x c
#'
dheve <- function(A, vertical) {
  if (vertical == TRUE) {
    elements <- A[rep(1:nrow(A), each=ncol(A)), ]
  } else {
    elements <- matrix(c(t(A)))[, rep(1, ncol(A))]
  }
  return(elements)
}


#' Aggregates elements of DHE and DVE matrices in a step to build
#' evidence matrix E.
#'
#' @param dhe DHE matrix of size Nc x c.
#' 
#' @param dve DVE matrix of size Nc x c.
#'
#' @return Matrix of size Nc x 1.
#'
gamma_fcm <- function(dhe, dve) {
  1 / rowSums(dhe/dve)
}


#' Rearranges elements of input matrix from a block matrix with vertical blocks
#' (column vectors) to a block matrix with horizontal blocks (row vectors).
#'
#' @param A Matrix of size Nc x 1.
#' 
#' @param c Number of columns in the wanted matrix.
#' Associated with the number of clusters.
#'
#' @return Matrix of size N x c.
#'
xi_fcm <- function(A, c) {
  matrix(A, ncol = c, byrow = TRUE)
}


#' Calculates data evidence matrix E from distances matrix D.
#'
#' @param D Distances matrix of size N x c.
#'
#' @return Matrix of size N x c.
#'
calculate_evidence <- function(D) {
  dve <- dheve(D, vertical = TRUE)
  dhe <- dheve(D, vertical = FALSE)
  output <- xi_fcm(A = gamma_fcm(dhe, dve), c = ncol(dhe))
  return(output)
}


#' Equation to calculate clusters' prototypes matrix \eqn{\hat{V}}.
#'
#' @param Phi Matrix with weights of size N x c.
#'
#' @param X Matrix with predictors of size N x p.
#'
#' @return Clusters' prototypes matrix of size c x p.
#'
estimate_V <- function(Phi, X) {
  V <- (t(Phi) %*% X) / colSums(Phi)
  return(V)
}


#' Estimated U matrix with memberships in semi-supervised case.
#'
#' @param D Distances matrix of size N x c.
#'
#' @param superF
#' Binary supervision matrix of size N x c.
#'
#' @param alpha
#' Scaling factor, a floating point > 0 regulating the impact of partial supervision.
#'
estimate_U <-
  function(
    D,
    superF,
    alpha
  ) {
    E <- calculate_evidence(D)
    
    U <- (E + alpha * superF) / (1 + alpha * rowSums(superF))
    
    return(U)
  }


#' Estimated T matrix with typicalities in unsupervised case.
#'
#' @param D Distances matrix of size N x c.
#' 
#' @param gammas
#' a c-vector of cluster-specific gamma hyperparameters.
#'
estimate_T <-
  function(
    D,
    gammas
  ) {
    G <- matrix(gammas, nrow = 1)[rep(1, nrow(D)), ]
    Tp <- G / (G + D)
    
    return(Tp)
  }


#' Initialization procedure to calculate values of gamma hyperparameters.
#' 
#' @param .model estimated model of class `fcm`
#' 
#' @param .X features matrix of size N x c
#'
init_gamma <- function(.model, .X) {
  # recreate distances
  .D <- .model$function_dist(.X, .model$V)
  out <- colSums(.model$U * .D) / colSums(.model$U)
  
  return(out)
}


#' Estimated T matrix with typicalities in semi-supervised case.
#'
#' @param D Distances matrix of size N x c.
#'
#' @param superF
#' Binary supervision matrix of size N x c.
#'
#' @param alpha
#' Scaling factor, a floating point > 0 regulating the impact of partial supervision.
#' 
#' @param gammas
#' a c-vector of cluster-specific gamma hyperparameters.
#' 
#' @param b
#' a scalar weighting the contribution of possibilistic membership in
#' SPFCM (semi-supervised possibilistic fuzzy c-means) model.
#' It is set to 1 by default for other semi-supervised models.
#' 
estimate_super_T <-
  function(
    D,
    superF,
    alpha,
    gammas,
    b = 1
  ) {
    G <- matrix(gammas, nrow = 1)[rep(1, nrow(D)), ]
    Tp <- (G + D * alpha * superF) / (G + D * (b + alpha * rowSums(superF)))
    
    return(Tp)
  }
