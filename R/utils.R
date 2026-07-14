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
    elements <- A[rep(1:nrow(A), each = ncol(A)), ]
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
  1 / rowSums(dhe / dve)
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


#' Adaptive Mahalanobis Distance for SPFCM
#'
#' Computes adaptive Mahalanobis distances between observations and cluster
#' prototypes using cluster-specific covariance matrices estimated from a
#' composite weight matrix.
#'
#' @param X Numeric features matrix of observations with dimensions \eqn{N \times p}.
#'
#' @param V Numeric matrix of cluster prototypes with dimensions
#' \eqn{c \times p}, where each row is the center of one cluster.
#'
#' @param Phi Numeric matrix of weights with dimensions
#' \eqn{N \times c}. Column \code{k} contains the observation weights used
#' to estimate the covariance matrix of cluster \code{k}.
#'
#' @param rho Optional numeric vector of length \eqn{c} containing the
#' determinant constraints for each cluster. Default is \code{rep(1, c)}.
#'
#' @return A numeric distances matrix of dimensions \eqn{N \times c}, where element
#' \code{D[j,k]} is the adaptive Mahalanobis distance between observation
#' \code{j} and prototype \code{k}.
#'
#' @details
#' The weight matrix \code{Phi} (\eqn{\Phi = [\phi_{jk}]}) corresponds to the observation-cluster weights
#' used in the adaptive metric estimation:
#'
#' \deqn{
#' \phi_{jk} =
#' a u_{jk}^{m}
#' +
#' b t_{jk}^{2}
#' +
#' \alpha b_j(t_{jk}-f_{jk})^{2}
#' }
#'
#' For each cluster, a weighted covariance matrix is computed:
#'
#' \deqn{
#' C_k =
#' \sum_j \phi_{jk}
#' (x_j-v_k)(x_j-v_k)^T
#' }
#'
#' The adaptive Mahalanobis distance is then:
#'
#' \deqn{
#' D_{jk} =
#' \rho_k^{1/p}
#' (\det(C_k))^{1/p}
#' (x_j-v_k)^T C_k^{-1}(x_j-v_k)
#' }
#'
#' where \code{p} is the number of features.
#'
#' The implementation avoids explicit matrix inversion and computes the
#' determinant scaling factor using a logarithmic determinant for numerical
#' stability.
#'
#' @references
#' Antoine et al. (2022) <doi:10.1016/j.fss.2022.08.003>
#' 
#' @export
#' 
spfcm_mahalanobis <- function(X, V, Phi, rho = NULL) {
  
  N <- nrow(X)
  p <- ncol(X)
  c <- nrow(V)
  
  if (is.null(rho)) {
    rho <- rep(1, c)
  }
  
  stopifnot(length(rho) == c)
  stopifnot(all(dim(Phi) == c(N, c)))
  
  D <- matrix(0, nrow = N, ncol = c)
  
  for (k in seq_len(c)) {
    
    aux <- sweep(X, 2, V[k, ], FUN = "-")
    w <- Phi[, k]
    
    Ck <- crossprod(aux, aux * w)
    
    DMah <- stats::mahalanobis(
      X,
      center = V[k, ],
      cov = Ck
    )
    
    log_det <- determinant(Ck, logarithm = TRUE)$modulus
    
    det_scale <- exp(as.numeric(log_det) / p)
    
    D[, k] <- rho[k]^(1 / p) * det_scale * DMah
  }
  
  return(D)
}
