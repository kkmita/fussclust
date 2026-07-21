set.seed(42)

X <- matrix(rnorm(100), ncol = 2)
superF <- matrix(0, nrow = nrow(X), ncol = ncol(X))
superF[1:10, 1] <- 1
superF[11:20, 2] <- 1

X_bigger <- rbind(X, X)

Xdummy <- matrix(c(rep(1, 20), rep(2, 20)), ncol = 2, byrow = TRUE)

Udummy <- matrix(stats::runif(40), ncol = 2)

U0 <- matrix(runif(nrow(X) * 3), ncol = 3)
