test_that("Number of clusters in SSFCM", {
  expect_error(SSFCM(X = X, C = 3, superF = superF, alpha = 1),
  "number of columns in `X` must match `C`.")
})

test_that("Correct alpha in SSFCM", {
  expect_error(SSFCM(X = X, C = 2, superF = superF, alpha = matrix(1)),
               "alpha must be either NULL or a scalar.")
})


test_that("Dimension of superF in SSFCM", {
  expect_error(SSFCM(X = X_bigger, C = 2, superF = superF, alpha = 1),
               "dimension of `superF` must be the same as dimension of `U`.")
})
