test_that("Correct gammas in PCM", {
  expect_error(
    PCM(X = X, C = 2, gammas = c(1, 1, 1)),
    "gammas must be either NULL or a numeric vector of length C."
  )
})
