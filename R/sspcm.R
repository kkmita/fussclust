#' Possibilistic c-Means model.
#'
#' @description
#' Unsupervised Possibilistic c-Means algorithm.
#' 
#' @param X
#' Features matrix *X*.
#'
#' @param C
#' Number of clusters.
#' 
#' @param gammas
#' Optionally: a vector of cluster-specific gamma hyperparameters.
#' Default value: `NULL`. In such case, the initialization value depends on
#' `initFCM` value. If `initFCM` is `NULL`, then vector filled with ones is used.
#' If `initFCM` is not `NULL`, then Fuzzy c-Means model is fitted, and the algorithm
#' implemented in `init_gamma` function is used to calculate the cluster-specific
#' gamma hyperparameters.
#'
#' @param initFCM
#' Default value: `NULL`. If `gammas` is not `NULL`, and `initFCM` is not `NULL`,
#' then special algorithm based on creating Fuzzy c-Means model and appropriately
#' calculating values of hyperparameters gamma is initiated.
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
#' @return An object of class `sspcm` containing:
#' \describe{
#'   \item{U}{An \eqn{N \times c} matrix of cluster memberships.}
#'   \item{V}{A \eqn{c \times p} matrix of cluster prototypes.}
#'   \item{function_dist}{An object of class `function` used to calculate distances.}
#'   \item{counter}{Integer number of iterations until convergence.}
#'   \item{gammas}{A vector of cluster-specific gamma hyperparameters.}
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
#' model_sspcm <- SSPCM(X = X, C = 2, superF = superF, alpha = 1)
#' print(model_sspcm$V)
#' 
SSPCM <- function(
    X,
    C,
    U = NULL,
    gammas = NULL,
    initFCM = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist,
    alpha = NULL,
    superF = NULL
) {
  if ( (!is.numeric(gammas) & !is.null(gammas)) 
        || (is.numeric(gammas) & length(gammas) != C)
     ){
    stop("gammas must be either NULL or a numeric vector of length C.", 
         call. = FALSE)
  }
  
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
  
  if (is.null(gammas)) {
    if (is.null(initFCM)) {
      gammas <- rep(1, C)
    } else {
      .modelFCM <- fussclust::FCM(X = X, C = C, U = U)
      gammas <- init_gamma(.modelFCM, X)
    }
  }
  
  counter = 0
  T_history <- list()
  V_history <- list()
  Phi_history <- list()
  
  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U
    
    Phi <- U_previous_iter^2 + (U_previous_iter-F)^2 * alpha * rowSums(superF)
    
    Phi_history[[counter]] <- Phi
    
    V <- estimate_V(Phi, X)
    
    V_history[[counter]] <- V
    
    U <- estimate_super_T(
      D = function_dist(X, V)^2,
      superF = superF,
      alpha = alpha,
      gammas = gammas
    )
    
    T_history[[counter]] <- U
    
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
    gammas = gammas,
    V_history = V_history,
    T_history = T_history,
    Phi_history = Phi_history,
    alpha = alpha
  )
  
  class(z) <- "sspcm"
  
  return(z)
}


#' Predict method for SSPCM objects
#'
#' @param object An object of class \code{sspcm}
#' 
#' @param X New data matrix of size \eqn{N \times p}
#' 
#' @param ... Not used
#'
#' @method predict sspcm
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
#' model_sspcm <- SSPCM(X = X, C = 2, superF = superF, alpha = 1)
#' predict(model_sspcm, matrix(rnorm(2), ncol = 2))
#' 
predict.sspcm <- function(object, X, ...) {
  Dpred <- object$function_dist(
    X,
    object$V
  )^2

  Tpred <- estimate_T(Dpred, object$gammas)

  return(Dpred)  
}