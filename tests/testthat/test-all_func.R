test_that("mylm works", {
  expect_equal(mylm(formula = Petal.Length ~ Sepal.Length + Sepal.Width, data = iris)$y[1], 1.4)
})


test_that("mylm works2", {
  expect_equal(mylm(formula = Petal.Length ~ Sepal.Length + Sepal.Width, data = iris)$X[1,2], 5.1)
})
