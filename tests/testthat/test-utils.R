
test_that("Distances Horizontally Exploded works", {
  expect_equal(dheve(matrix(c(5, 25, 61, 145, 85, 41), ncol = 2), vertical = FALSE),
               matrix(rep(c(5, 145, 25, 85, 61, 41), 2), ncol = 2))
})

test_that("Distances Vertically Exploded works", {
  expect_equal(dheve(matrix(c(5, 25, 61, 145, 85, 41), ncol = 2), vertical = TRUE),
               matrix(rep(c(5, 25, 61, 145, 85, 41), each = 2), ncol = 2))
})

