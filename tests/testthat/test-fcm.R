test_that("Prototypes matrix in dummy case in FCM", {
  expect_equal(
    FCM(X = Xdummy, U = Udummy, C = 2)$V,
    matrix(c(2, 2, 1, 1), ncol = 2, byrow = TRUE)
  )
})
