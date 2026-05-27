set.seed(42)

X <- matrix(rnorm(100), ncol = 2)
superF <- matrix(0, nrow = nrow(X), ncol = ncol(X))
superF[1:10, 1] <- 1
superF[11:20, 2] <- 1

X_bigger <- rbind(X, X)