test_that("Correct gammas in PCM", {
  expect_error(
    PCM(X = X, C = 2, gammas = c(1, 1, 1)),
    "gammas must be either NULL or a numeric vector of length C."
  )
})

test_that("No history stored in PCM", {
  expect_equal(
    PCM(X = X, C = 2)$U_history,
    NULL
  )
})

test_that("History store in PCM", {
  expect_equal(
    class(PCM(X = X, C = 2, store_history = TRUE)$U_history),
    "list"
  )
})

test_that("Different value of m changes results", {
  expect_false(
    isTRUE(
      all.equal(
        PCM(X = X, U = U0, C = 3, m = 2)$V,
        PCM(X = X, U = U0, C = 3, m = 10)$V
      )
    )
  )
})

test_that("Value of m must be > 1", {
  expect_error(
    PCM(X = X, C = 2, m = 1)
  )
})
