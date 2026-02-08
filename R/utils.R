#' Creates DHE (stands for "distances horizontally exploded") and DVE
#' (stands for "distances vertically exploded") matrices.
#'
#' @param A Matrix of size N x c.
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
#' @param dve DVE matrix of size Nc x c.
#'
#' @return Matrix of size Nc x 1.
#'
gamma_fcm <- function(dhe, dve) {
  # 1 / ((dhe/dve) %*% matrix(rep(1, ncol(dhe))))
  1 / rowSums(dhe/dve)
}


#' Rearranges elements of input matrix from a block matrix with vertical blocks
#' (column vectors) to a block matrix with horizontal blocks (row vectors).
#'
#' @param A Matrix of size Nc x 1.
#' @param c Number of columns in the wanted matrix.
#' Associated with the number of clusters.
#'
#' @return Matrix of size N x c.
#'
xi_fcm <- function(A, c) {
  matrix(A, ncol=c, byrow=TRUE)
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


#' Equation to calculate clusters' prototypes matrix $\hat{V}$.
#'
#' @param Phi Matrix with weights of size N x c.
#'
#' @param X Matrix with predictors of size N x p.
#'
#' @return Clusters' prototypes matrix of size c x p.
#' @export
#'
estimate_V <- function(Phi, X) {
  V <- (t(Phi) %*% X) / colSums(Phi)
  return(V)
}


#' Estimated U matrix with memberships.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (c, p)
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
estimate_U <-
  function(
    X,
    V,
    superF,
    alpha,
    function_dist
  ) {
    D <- function_dist(X, V)^2
    E <- calculate_evidence(D)
    
    U <- (E + alpha * superF) / (1 + alpha * rowSums(superF))
    
    return(U)
  }



#' Estimated T matrix with typicalities.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (c, p)
#'
#' @param superF
#' the supervisionbinary matrix of the same dimension as *U*.
#'
#' @param alpha
#' a *N*-vector of observation-specific scaling factor values
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param gammas
#' a *c*-vector of cluster-specific gamma parameter
#'
estimate_T <-
  function(
    X,
    V,
    function_dist,
    gammas
  ) {
    D <- function_dist(X, V)^2
    G <- matrix(gammas, nrow = 1)[rep(1, nrow(D)), ]
    Tp <- G / (G + D)
    
    return(Tp)
  }


init_gamma <- function(.model, .X) {
  # recreate distances
  .D <- .model$function_dist(.X, .model$V)
  out <- colSums(.model$U * .D) / colSums(.model$U)
  
  return(out)
}


#' Estimated T matrix with typicalities.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (c, p)
#'
#' @param superF
#' the supervisionbinary matrix of the same dimension as *U*.
#'
#' @param alpha
#' a *N*-vector of observation-specific scaling factor values
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param gammas
#' a *c*-vector of cluster-specific gamma parameter
#'
estimate_super_T <-
  function(
    X,
    V,
    superF,
    alpha,
    function_dist,
    gammas
  ) {
    D <- function_dist(X, V)^2
    G <- matrix(gammas, nrow = 1)[rep(1, nrow(D)), ]
    Tp <- (G + D * alpha * superF) / (G + D * (1 + alpha * rowSums(superF)))
    
    return(Tp)
  }
