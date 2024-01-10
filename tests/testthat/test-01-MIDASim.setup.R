# context("Testing MIDASim.setup()")

test_that("MIDASim.setup works properly for setting up MIDASim nonparametric model", {

  data(count.ibd)
  fitted.nonparametric = MIDASim.setup(count.ibd, mode = "nonparametric")

  expect_equal(dim(fitted.nonparametric$mat01), dim(count.ibd))
  expect_true("rel.abund.1" %in% names(fitted.nonparametric))

})

test_that("MIDASim.setup works properly for setting up MIDASim parametric model", {

  data(count.ibd)
  fitted.parametric = MIDASim.setup(count.ibd, mode = "parametric")

  expect_equal(dim(fitted.parametric$mat01), dim(count.ibd))
  expect_equal(length(fitted.parametric$mu.est), ncol(count.ibd))
  expect_true("mu.est" %in% names(fitted.parametric))
  expect_true("sigma.est" %in% names(fitted.parametric))
  expect_true("Q.est" %in% names(fitted.parametric))

})
