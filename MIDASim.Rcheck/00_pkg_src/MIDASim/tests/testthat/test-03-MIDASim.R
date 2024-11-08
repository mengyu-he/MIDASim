# context("Testing MIDASim()")

test_that("MIDASim simulates data with correct dimension", {

  data(count.ibd)
  fitted.nonparametric = MIDASim.setup(count.ibd, mode = "nonparametric")
  fitted.parametric = MIDASim.setup(count.ibd, mode = "parametric")

  fitted.nonparametric.modified = MIDASim.modify(fitted.nonparametric)
  fitted.parametric.modified = MIDASim.modify(fitted.parametric)

  sim.nonparametric = MIDASim(fitted.nonparametric.modified)
  sim.parametric = MIDASim(fitted.parametric.modified)

  # nonparametric mode
  expect_equal(dim(sim.nonparametric$sim_count), dim(count.ibd))

  # parametric mode
  expect_equal(dim(sim.parametric$sim_count), dim(count.ibd))

})
