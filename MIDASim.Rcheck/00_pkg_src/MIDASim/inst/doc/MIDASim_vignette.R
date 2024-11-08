## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  eval = FALSE)

## ----load package, echo=FALSE-------------------------------------------------
#  library(MIDASim)
#  library(MASS)
#  library(psych, quietly = TRUE)
#  library(scam, quietly = TRUE)
#  library(pracma, quietly = TRUE)
#  library(vegan, quietly = TRUE)

## -----------------------------------------------------------------------------
#  devtools::install_local("MIDASim_0.1.0.tar.gz", dependencies = TRUE)

## -----------------------------------------------------------------------------
#  install_github("mengyu-he/MIDASim")

## -----------------------------------------------------------------------------
#  browseVignettes("MIDASim")

## -----------------------------------------------------------------------------
#  vignette("MIDASim_vignette", package = "MIDASim")

## -----------------------------------------------------------------------------
#  data("count.ibd")

## ----eval = FALSE-------------------------------------------------------------
#  library(HMP2Data)
#  IBD16S()
#  
#  count.ibd = t(IBD16S_mtx[, colSums(IBD16S_mtx)>3000] )
#  count.ibd = count.ibd[, colSums(count.ibd>0)>1]
#  
#  dim(count.ibd)

## -----------------------------------------------------------------------------
#  count.ibd.setup = MIDASim.setup(otu.tab, mode = 'nonparametric', n.break.ties = 10)

## -----------------------------------------------------------------------------
#  count.ibd.setup = MIDASim.setup(otu.tab, mode = 'parametric')

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                     lib.size = NULL,
#                                     mean.rel.abund = NULL,
#                                     gengamma.mu = NULL,
#                                     sample.1.prop = NULL,
#                                     taxa.1.prop = NULL)

## -----------------------------------------------------------------------------
#  new.lib.size = sample(1000:10000, size=700, replace=TRUE)
#  
#  obs.sample.1.ct = rowSums(count.ibd > 0)
#  xvar = log10(rowSums(count.ibd))
#  scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
#  sample.1.ct = 10^(predict(scamfit.non0,
#                             newdata = data.frame(xvar = log10(new.lib.size) )) )
#  
#  n.taxa = ncol(count.ibd)
#  new.sample.1.prop = sample.1.ct/n.taxa
#  new.taxa.1.prop = fitted$taxa.1.prop * (sum(new.sample.1.prop) * n.taxa / 700) / sum(fitted$taxa.1.prop)

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      lib.size = new.lib.size,
#                                      sample.1.prop = new.sample.1.prop,
#                                      taxa.1.prop = new.taxa.1.prop)

## -----------------------------------------------------------------------------
#  beta = 0.1
#  new.mean.rel.abund = count.ibd.setup$mean.rel.abund
#  new.mean.rel.abund[1:10] = exp(beta) * new.mean.rel.abund[1:10]
#  new.mean.rel.abund = new.mean.rel.abund / sum(new.mean.rel.abund)

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      mean.rel.abund = new.mean.rel.abund)

## -----------------------------------------------------------------------------
#  new.taxa.1.prop = count.ibd.setup$taxa.1.prop
#  new.taxa.1.prop[1:10] = new.taxa.1.prop[1:10] + 0.05

## -----------------------------------------------------------------------------
#  r = sum(new.taxa.1.prop) / sum(count.ibd.setup$taxa.1.prop)
#  new.sample.1.prop = count.ibd.setup$sample.1.prop * r

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      taxa.1.prop = new.taxa.1.prop,
#                                      sample.1.prop = new.sample.1.prop)

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      mean.rel.abund = new.mean.rel.abund,
#                                      taxa.1.prop = new.taxa.1.prop,
#                                      sample.1.prop = new.sample.1.prop)

## -----------------------------------------------------------------------------
#  new.lib.size = sample(1000:10000, size=700, replace=TRUE)
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      lib.size = new.lib.size)

## -----------------------------------------------------------------------------
#  beta = 0.1
#  new.mean.rel.abund = count.ibd.setup$mean.rel.abund
#  new.mean.rel.abund[1:10] = exp(beta) * new.mean.rel.abund[1:10]
#  new.mean.rel.abund = new.mean.rel.abund / sum(new.mean.rel.abund)

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      mean.rel.abund = new.mean.rel.abund)

## -----------------------------------------------------------------------------
#  count.ibd.modified = MIDASim.modify(count.ibd.setup,
#                                      gengamma.mu = count.ibd.setup$mu.est + beta)

## -----------------------------------------------------------------------------
#  simulated.data = MIDASim(count.ibd.modified)
#  summary(simulated.data)

