#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' Features matrix *X*.
#'
#' @param C
#' Number of clusters.
#'
#' @param U
#' Optionally: a concrete initialization memberships matrix.
#' Used mainly for reproducibility.
#' Default value `NULL` - algorithm uses random initialization in such case.
#'
#' @param max_iter
#' Maximum number of iterations. Default value: 200.
#' 
#' @param conv_criterion
#' Convergence criterion value used at the end of each iteration of
#' Alternating Optimization algorithm.
#'
#' @param function_dist
#' Optionally: a function of two arguments: matrices *X* and *V* of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of *X* and all rows of *V*.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' Scaling factor, a floating point > 0 regulating the impact of partial supervision.
#'
#' @param superF
#' Binary supervision matrix of the same dimension as *U*.
#'
#' @export
#' 
#' @return An object of class `ssfcm` containing:
#' \describe{
#'   \item{U}{An \eqn{N \times c} matrix of cluster memberships.}
#'   \item{V}{A \eqn{c \times p} matrix of cluster prototypes.}
#'   \item{function_dist}{An object of class `function` used to calculate distances.}
#'   \item{counter}{Integer number of iterations until convergence.}
#'   \item{V_history}{A list of length `counter` with \eqn{c \times p} 
#'   prototypes matrices estimated in each loop of the algorithm.}
#'   \item{U_history}{A list of length `counter` with \eqn{N \times c} 
#'   memberships matrices estimated in each loop of the algorithm.}
#'   \item{Phi_history}{A list of length `counter` with \eqn{N \times c} 
#'   phi weights in each loop of the algorithm.}
#'   \item{alpha}{A value of the scaling factor used.}
#' }
#' 
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#' superF <- matrix(0, nrow = nrow(X), ncol = ncol(X))
#' superF[1:10, 1] <- 1
#' superF[11:20, 2] <- 1
#' model_ssfcm <- SSFCM(X = X, C = 2, superF = superF, alpha = 1)
#' print(model_ssfcm$V)
#' 
SSFCM <- function(
    X,
    C,
    U = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist,
    alpha = NULL,
    superF = NULL
) {
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.null(dim(alpha))) {
    stop("alpha must be either NULL or a scalar.", call. = FALSE)
  }
  
  if (ncol(X) != C) {
    stop("number of columns in `X` must match `C`.", call. = FALSE)
  }
  
  if (!is.null(superF) && ( (nrow(superF) != nrow(X)) || (ncol(superF) != C) )) {
    stop("dimension of `superF` must be the same as dimension of `U`.")
  }  
  
  if (is.null(U)) {
    U <- matrix(stats::runif(nrow(X)*C), ncol=C)
  }
  
  # Rows of U should sum up to 1
  U <- t(apply(U, 1, function(x) x / sum(x)))
  
  counter = 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()
  
  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U
    
    Phi <- U_previous_iter^2 + (U_previous_iter-superF)^2 * alpha * rowSums(superF)
    
    Phi_history[[counter]] <- Phi
    
    V <- estimate_V(Phi, X)
    
    V_history[[counter]] <- V
    
    U <- estimate_U(
      D = function_dist(X, V)^2,
      superF = superF,
      alpha = alpha)
    
    U_history[[counter]] <- U
    
    conv_iter <- base::norm(U - U_previous_iter, type="F")
    
    if (conv_iter < conv_criterion) {
      break
    }
  }
  
  z <- list(
    U = U,
    V = V,
    function_dist = function_dist,
    counter = counter,
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history,
    alpha = alpha
  )
  
  class(z) <- "ssfcm"
  
  return(z)
}


#' Predict method for SSFCM objects
#'
#' @param object An object of class \code{ssfcm}
#' 
#' @param X New data matrix of size \eqn{N \times p}
#' 
#' @param ... Not used
#'
#' @method predict ssfcm
#' 
#' @export
#' 
#' @return A matrix of size \eqn{N \times C}, where \eqn{C} is the number of columns
#' in `object$U` containing predicted memberships.
#' 
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#' superF <- matrix(0, nrow = nrow(X), ncol = ncol(X))
#' superF[1:10, 1] <- 1
#' superF[11:20, 2] <- 1
#' model_ssfcm <- SSFCM(X = X, C = 2, superF = superF, alpha = 1)
#' predict(model_ssfcm, matrix(rnorm(2), ncol = 2))
predict.ssfcm <- function(object, X, ...) {
  Dpred <- object$function_dist(
    X,
    object$V
  )^2
  
  Upred <- calculate_evidence(Dpred)
  
  return(Upred)
}
