#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' a matrix *X* with predictor variables.
#'
#' @param C
#' a number of clusters to find.
#'
#' @param U
#' optionally: a first memberships matrix to initialize the algorithm.
#' Used mainly for reproducibility to compare calculations with other packages
#' (e.g. in Python).
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @export
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
  
  if (!is.null(alpha) && length(alpha) != 1) {
    stop("'alpha' must be either NULL or a scalar (length 1 value).")
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
      X = X,
      V = V,
      superF = superF,
      alpha = alpha,
      function_dist = function_dist)
    
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
#' @param X New data matrix
#' @param ... Not used
#'
#' @return A data frame with predicted labels
#'
#' @method predict ssfcm
#' @export
predict.ssfcm <- function(object, .X, ...) {
  Dpred <- object$function_dist(
    .X,
    object$V
  )^2
  
  Upred <- calculate_evidence(Dpred)
  
  return(Upred)
}
