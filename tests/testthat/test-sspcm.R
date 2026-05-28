test_that("Number of clusters in SSPCM", {
  expect_error(
    SSPCM(X = X, C = 3, superF = superF, alpha = 1),
    "number of columns in `superF` must match `C`."
  )
})

test_that("Correct alpha in SSPCM", {
  expect_error(
    SSPCM(X = X, C = 2, superF = superF, alpha = matrix(1)),
    "alpha must be either NULL or a scalar."
  )
})

test_that("Dimension of superF in SSPCM", {
  expect_error(
    SSPCM(X = X_bigger, C = 2, superF = superF, alpha = 1),
    "dimension of `superF` must be the same as dimension of `U`."
  )
})

test_that("Correct gammas in SSPCM", {
  expect_error(
    SSPCM(X = X, C = 2, superF = superF, alpha = 1, gammas = c(1, 1, 1)),
    "gammas must be either NULL or a numeric vector of length C."
  )
})

test_that("Correct prediction format in SSPCM", {
  expect_equal(
    dim(predict(SSPCM(X = X, C = 2, superF = superF, alpha = 1), matrix(rnorm(4), ncol = 2))),
    c(2, 2)
  )
})
