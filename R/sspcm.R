#' @export
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
    
    # Phi <- Tm_previous_iter^2
    
    Phi <- U_previous_iter^2 + (U_previous_iter-F)^2 * alpha * rowSums(superF)
    
    
    # Modify `Phi` if running semi-supervised PCM
    # if (!is.null(alpha)) {
    #    Tm_alpha <- alpha * (Tm_previous_iter - superF)^2
    #   Phi <- Phi + Tm_alpha
    # }
    
    Phi_history[[counter]] <- Phi
    
    V <- estimate_V(Phi, X)
    
    V_history[[counter]] <- V
    
    U <- estimate_super_T(
      X = X,
      V = V,
      superF = superF,
      alpha = alpha,
      function_dist = function_dist,
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
    Phi_history = Phi_history
  )
  
  class(z) <- "sspcm"
  
  return(z)
}



#' Soft assignment score function
#'
#' @param object object
#' @param newdata data to be predicted
#'
#' @return typicalities matrix
#'
#' @export
predict.sspcm <- function(object, newdata) {
  output <- estimate_T(
    X = newdata,
    V = object$V,
    superF = NULL,
    alpha = NULL,
    function_dist = object$function_dist,
    gammas = object$gammas
  )
  return(output)
}