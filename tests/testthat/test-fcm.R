test_that("Prototypes matrix in dummy case in FCM", {
  expect_equal(
    FCM(X = Xdummy, U = Udummy, C = 2)$V,
    matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE)
  )
})

test_that("No history stored in FCM", {
  expect_equal(
    FCM(X = X, C = 2)$U_history,
    NULL
  )
})

test_that("History store in FCM", {
  expect_equal(
    class(FCM(X = X, C = 2, store_history = TRUE)$U_history),
    "list"
  )
})
