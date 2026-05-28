#' Fuzzy C-Means clustering model
#'
#' @description
#' Fits a Fuzzy C-Means (FCM) clustering model using the Alternating
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
#' @return An object of class `fcm` containing:
#' \describe{
#'   \item{U}{An \eqn{N \times C} membership matrix.}
#'   \item{V}{A \eqn{C \times p} matrix of cluster prototypes.}
#'   \item{function_dist}{The distance function used by the model.}
#'   \item{counter}{Number of iterations performed until convergence.}
#'   \item{V_history}{A list of length `counter` containing prototype
#'   matrices estimated at each iteration.}
#'   \item{U_history}{A list of length `counter` containing membership
#'   matrices estimated at each iteration.}
#'   \item{Phi_history}{A list of length `counter` containing phi-weight
#'   matrices estimated at each iteration.}
#' }
#' 
#' @references 
#' Bezdek, J. C. (1981). 
#' \emph{Pattern Recognition with Fuzzy Objective Function Algorithms}. 
#' Springer US.
#' https://doi.org/10.1007/978-1-4757-0450-1
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#'
#' model_fcm <- fussclust::FCM(
#'   X = X,
#'   C = 2
#' )
#'
#' print(model_fcm$V)
#'
#' @export
FCM <- function(
  X,
  C,
  U = NULL,
  max_iter = 200,
  conv_criterion = 1e-4,
  function_dist = rdist::cdist
) {
  if (is.null(U)) {
    U <- matrix(stats::runif(nrow(X) * C), ncol = C)
  }

  # Rows of U should sum up to 1
  U <- t(apply(U, 1, function(x) x / sum(x)))

  counter <- 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()

  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    Phi <- U_previous_iter^2

    Phi_history[[counter]] <- Phi

    V <- estimate_V(Phi, X)

    V_history[[counter]] <- V

    D <- function_dist(X, V)^2
    U <- calculate_evidence(D)

    U_history[[counter]] <- U

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
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history
  )

  class(z) <- "fcm"

  return(z)
}
