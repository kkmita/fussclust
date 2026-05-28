#' Possibilistic C-Means clustering model
#'
#' @description
#' Fits a Possibilistic C-Means (PCM) clustering model using the Alternating
#' Optimization algorithm.
#'
#' @param X A numeric feature matrix.
#'
#' @param C Integer specifying the number of clusters.
#'
#' @param U Optional initial membership matrix.
#' Primarily intended for reproducibility purposes.
#' If `NULL` (default), the algorithm uses a random initialization.
#'
#' @param gammas Optional vector of cluster-specific gamma hyperparameters.
#' If `NULL` (default), the initialization strategy depends on the value
#' of `initFCM`.
#'
#' If `initFCM` is `NULL`, a vector of ones is used.
#' Otherwise, a Fuzzy C-Means model is first fitted, and the
#' `init_gamma()` function is used to estimate the cluster-specific
#' gamma hyperparameters.
#'
#' @param initFCM Optional fitted Fuzzy C-Means model used to initialize
#' cluster-specific gamma hyperparameters via weighted averaging.
#' If `NULL` (default), no preliminary Fuzzy C-Means initialization is used.
#' If provided, this argument is effective only when `gammas` is `NULL`.
#'
#' @param max_iter Maximum number of iterations.
#' Defaults to `200`.
#'
#' @param conv_criterion Convergence threshold used at the end of each
#' iteration of the Alternating Optimization algorithm.
#'
#' @param function_dist Optional distance function.
#' The function must accept two matrices, `X` and `V`, with the same
#' number of columns, and return a matrix of size
#' `nrow(X) x nrow(V)` containing distances between each row of `X`
#' and each row of `V`.
#'
#' For the Euclidean distance, the returned distances should not be squared.
#' Defaults to [rdist::cdist()].
#'
#' @param store_history Logical indicating whether optimization
#' histories should be stored. If `FALSE`, the returned object
#' will contain `NULL` history fields. Defaults to `TRUE`.
#'
#' @return An object of class `pcm` containing:
#' \describe{
#'   \item{U}{An \eqn{N \times C} membership matrix.}
#'   \item{V}{A \eqn{C \times p} matrix of cluster prototypes.}
#'   \item{function_dist}{The distance function used by the model.}
#'   \item{counter}{Number of iterations performed until convergence.}
#'   \item{gammas}{Vector of cluster-specific gamma hyperparameters.}
#'   \item{U_history}{If `store_history = TRUE`, a list of length
#'   `counter` containing membership matrices estimated at each
#'   iteration; otherwise `NULL`.}
#'   \item{V_history}{If `store_history = TRUE`, a list of length
#'   `counter` containing prototype matrices estimated at each
#'   iteration; otherwise `NULL`.}
#'   \item{Phi_history}{If `store_history = TRUE`, a list of length
#'   `counter` containing phi-weight matrices estimated at each
#'   iteration; otherwise `NULL`.}
#' }
#'
#' @references
#' Krishnapuram, R., & Keller, J. (1993).
#' \emph{A possibilistic approach to clustering}.
#' IEEE Transactions on Fuzzy Systems, 1(2), 98–110.
#' https://doi.org/10.1109/91.227387
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#'
#' model_pcm <- fussclust::PCM(
#'   X = X,
#'   C = 2,
#'   initFCM = TRUE
#' )
#'
#' print(model_pcm$V)
#'
#' @export
PCM <- function(
  X,
  C,
  U = NULL,
  gammas = NULL,
  initFCM = NULL,
  max_iter = 200,
  conv_criterion = 1e-4,
  function_dist = rdist::cdist,
  store_history = FALSE
) {
  if ((!is.numeric(gammas) & !is.null(gammas)) ||
    (is.numeric(gammas) & length(gammas) != C)
  ) {
    stop("gammas must be either NULL or a numeric vector of length C.",
      call. = FALSE
    )
  }

  if (is.null(U)) {
    U <- matrix(stats::runif(nrow(X) * C), ncol = C)
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

  counter <- 0

  if (store_history) {
    U_history <- list()
    V_history <- list()
    Phi_history <- list()
  } else {
    U_history <- NULL
    V_history <- NULL
    Phi_history <- NULL
  }

  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    Phi <- U_previous_iter^2
    V <- estimate_V(Phi, X)
    U <- estimate_T(D = function_dist(X, V)^2, gammas)

    if (store_history) {
      U_history[[counter]] <- U
      V_history[[counter]] <- V
      Phi_history[[counter]] <- Phi
    }

    conv_iter <- base::norm(U - U_previous_iter, type = "F")

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
    U_history = U_history,
    V_history = V_history,
    Phi_history = Phi_history
  )

  class(z) <- "pcm"

  return(z)
}
