test_that("Number of clusters in FCM", {
  expect_error(FCM(X = X, C = 3),
               "number of columns in `X` must match `C`.")
})
