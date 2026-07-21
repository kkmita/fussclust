#' Adaptive Mahalanobis Distance for Fuzzy Clustering
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
#' D. E. Gustafson and W. C. Kessel (1978) <doi: 10.1109/CDC.1978.268028>
#' Antoine et al. (2022) <doi:10.1016/j.fss.2022.08.003>
#'
#' @export
#'
adist <- function(X, V, Phi, rho = NULL) {
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
    aux <- sweep(X, 2, V[k,], FUN = "-")
    w <- Phi[, k]
    
    Ck <- crossprod(aux, aux * w)
    
    DMah <- stats::mahalanobis(X,
                               center = V[k,],
                               cov = Ck)
    
    log_det <- determinant(Ck, logarithm = TRUE)$modulus
    
    det_scale <- exp(as.numeric(log_det) / p)
    
    D[, k] <- rho[k]^(1 / p) * det_scale * DMah
  }
  
  return(D)
}
