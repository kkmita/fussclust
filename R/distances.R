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
#' D. E. Gustafson and W. C. Kessel (1978) <doi:10.1109/CDC.1978.268028>
#' 
#' @references
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


#' Control parameters for built-in distance functions
#'
#' @param rho Optional tuning parameter for adaptive distance.
#'
#' @return A list of control parameters.
#'
#' @export
#'
distance_control <- function(rho = NULL) {
  list(rho = rho)
}


#' Compute a distance matrix
#'
#' Computes the raw distance matrix using either one of the built-in
#' implementations or a user-defined distance function.
#'
#' Built-in implementations:
#'   - "cdist"
#'   - "adaptive"
#'
#' User-defined distance functions must accept a single argument `ctx`,
#' a named list containing the variables available in the current model.
#' 
#' Note that `rdist::cdist` returns \eqn{D = [d_{jk}]}, and hence its
#' result must be raised to the power of \eqn{2}, whereas `adist`
#' return as \eqn{D = [d^2_{jk}]}.
#'
#' @param ctx A named list representing the model context.
#' @param distance Either a built-in distance name or a custom function.
#' @param control Control parameters for built-in distance methods.
#'
#' @return A distance matrix.
#'
compute_distance <- function(ctx,
                             distance = "cdist",
                             control = distance_control()) {
  if (is.character(distance)) {
    distance <- match.arg(distance, c("cdist", "adaptive"))
    
    return(switch(
      distance,
      cdist = rdist::cdist(ctx$X, ctx$V, metric = "euclidean")^2,
      adaptive = adist(
        X = ctx$X,
        V = ctx$V,
        Phi = ctx$Phi,
        rho = control$rho
      )
    ))
  }
  
  if (is.function(distance)) {
    return(distance(ctx))
  }
  
  stop("`distance` must be one of 'cdist', 'adaptive', or a function.",
       call. = FALSE)
}