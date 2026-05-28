test_that("Correct gammas in PCM", {
  expect_error(
    PCM(X = X, C = 2, gammas = c(1, 1, 1)),
    "gammas must be either NULL or a numeric vector of length C."
  )
})

test_that("No history stored in PCM", {
  expect_equal(
    PCM(X = X, C = 2)$U_history,
    NULL)
})

test_that("History store in PCM", {
  expect_equal(
    class(PCM(X = X, C = 2, store_history = TRUE)$U_history),
    "list")
})