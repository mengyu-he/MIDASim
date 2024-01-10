# context("Testing MIDASim.modify()")

data(count.ibd)
fitted.nonparametric <- MIDASim.setup(count.ibd, mode = "nonparametric")
fitted.parametric <- MIDASim.setup(count.ibd, mode = "parametric")


test_that("eta and theta are correctly computed without modifications", {

  # nonparametric mode
  fitted.nonparametric.modified = MIDASim.modify(fitted.nonparametric)

  expect_equal(sum(is.finite(fitted.nonparametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.nonparametric.modified$theta)), ncol(count.ibd))

  # parametric mode
  fitted.parametric.modified = MIDASim.modify(fitted.parametric)

  expect_equal(sum(is.finite(fitted.parametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.parametric.modified$theta)), ncol(count.ibd))

})

test_that("eta and theta are correctly computed when lib.size is changed", {

  new.lib.size = sample(1000:10000, size=700, replace=T)

  # nonparametric mode
  obs.sample.1.ct = rowSums(count.ibd > 0)
  xvar = log10(rowSums(count.ibd))
  scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
  sample.1.ct = 10^(predict(scamfit.non0,
                            newdata = data.frame(xvar = log10(new.lib.size) )) )

  n.taxa = ncol(count.ibd)
  new.sample.1.prop = sample.1.ct/n.taxa
  new.taxa.1.prop = fitted.nonparametric$taxa.1.prop * (sum(new.sample.1.prop) * n.taxa / 700) / sum(fitted.nonparametric$taxa.1.prop)

  fitted.nonparametric.modified = MIDASim.modify(fitted.nonparametric,
                                                 lib.size = new.lib.size,
                                                 sample.1.prop = new.sample.1.prop,
                                                 taxa.1.prop = new.taxa.1.prop)

  expect_equal(sum(is.finite(fitted.nonparametric.modified$eta)), length(new.lib.size))
  expect_equal(sum(is.finite(fitted.nonparametric.modified$theta)), ncol(count.ibd))


  # parametric mode
  fitted.parametric.modified = MIDASim.modify(fitted.parametric,
                                              lib.size = new.lib.size)
  expect_equal(sum(is.finite(fitted.parametric.modified$eta)), length(new.lib.size))
  expect_equal(sum(is.finite(fitted.parametric.modified$theta)), ncol(count.ibd))

})

test_that("eta and theta are correctly computed when mean.rel.abund is changed", {

  beta = 0.1
  new.mean.rel.abund = fitted.nonparametric$mean.rel.abund
  new.mean.rel.abund[1:10] = exp(beta) * new.mean.rel.abund[1:10]
  new.mean.rel.abund = new.mean.rel.abund / sum(new.mean.rel.abund)

  # nonparametric mode
  fitted.nonparametric.modified = MIDASim.modify(fitted.nonparametric,
                                                 mean.rel.abund = new.mean.rel.abund)

  expect_equal(sum(is.finite(fitted.nonparametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.nonparametric.modified$theta)), ncol(count.ibd))

  # parametric mode
  fitted.parametric.modified = MIDASim.modify(fitted.parametric,
                                              mean.rel.abund = new.mean.rel.abund)

  expect_equal(sum(is.finite(fitted.parametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.parametric.modified$theta)), ncol(count.ibd))


})

test_that("In nonparametric mode, eta and theta are correctly computed when taxa.1.prop and sample.1.prop are changed", {

  new.taxa.1.prop = fitted.nonparametric$taxa.1.prop
  new.taxa.1.prop[1:10] = new.taxa.1.prop[1:10] + 0.05
  r = sum(new.taxa.1.prop) / sum(fitted.nonparametric$taxa.1.prop)
  new.sample.1.prop = fitted.nonparametric$sample.1.prop * r

  fitted.nonparametric.modified = MIDASim.modify(fitted.nonparametric,
                                                 taxa.1.prop = new.taxa.1.prop,
                                                 sample.1.prop = new.sample.1.prop)

  expect_equal(sum(is.finite(fitted.nonparametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.nonparametric.modified$theta)), ncol(count.ibd))

})

test_that("In parametric mode, eta and theta are correctly computed when gengamma.mu is changed", {

  beta = 0.1
  fitted.parametric.modified = MIDASim.modify(fitted.parametric,
                                              gengamma.mu = fitted.parametric$mu.est + beta)

  expect_equal(sum(is.finite(fitted.parametric.modified$eta)), nrow(count.ibd))
  expect_equal(sum(is.finite(fitted.parametric.modified$theta)), ncol(count.ibd))

})
