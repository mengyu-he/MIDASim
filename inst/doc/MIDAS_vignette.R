## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  eval = FALSE)

## ----load package, echo=FALSE-------------------------------------------------
#  library(MIDAS)
#  library(MASS)
#  library(psych, quietly = TRUE)
#  library(scam, quietly = TRUE)
#  library(pracma, quietly = TRUE)
#  library(vegan, quietly = TRUE)

## -----------------------------------------------------------------------------
#  devtools::install_local("MIDAS_0.1.0.tar.gz", dependencies = TRUE)

## -----------------------------------------------------------------------------
#  install_github("mengyu-he/MIDAS")

## -----------------------------------------------------------------------------
#  browseVignettes("MIDAS")

## -----------------------------------------------------------------------------
#  vignette("MIDAS_vignette", package = "MIDAS")

## -----------------------------------------------------------------------------
#  data("count.ibd")

## -----------------------------------------------------------------------------
#  library(HMP2Data)
#  IBD16S()

## -----------------------------------------------------------------------------
#  data("count.ibd")

## ---- eval = FALSE------------------------------------------------------------
#  library(HMP2Data)
#  IBD16S()
#  
#  count.ibd <- t(IBD16S_mtx[, colSums(IBD16S_mtx)>3000] )
#  count.ibd <- count.ibd[, colSums(count.ibd>0)>1]
#  
#  dim(count.ibd)

## -----------------------------------------------------------------------------
#  fitted <- Midas.setup(count.ibd, n.break.ties = 100, fit.beta = F)

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                 lib.size = NULL,
#                                 rel.abund = NULL,
#                                 taxa.1.prop = NULL)

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  lib.size = sample(1000:10000,
#                                                    nrow(count.ibd),
#                                                    replace = T) )

## -----------------------------------------------------------------------------
#  n.taxa <- ncol(count.ibd)
#  target.rel.abund <- NULL
#  
#  for (j in 1:n.taxa) {
#    target.rel.abund[j] <- -1
#    while (target.rel.abund[j]<0 | target.rel.abund[j]>1) {
#      target.rel.abund[j] <- fitted$obs.rel.abund[j] + 0.1*(rnorm(1, fitted$obs.rel.abund[j], fitted$obs.rel.abund[j]))
#    }
#  }

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  rel.abund  = target.rel.abund)

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  taxa.1.prop = sample(fitted$prop.1))

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  rel.abund = target.rel.abund,
#                                  taxa.1.prop = "same")

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  lib.size = sample(1000:10000, nrow(count.ibd), replace = T),
#                                  taxa.1.prop  = sample(fitted$prop.1),
#                                  method.total0 = "scam")

## -----------------------------------------------------------------------------
#  fitted.modified <- Midas.modify(fitted,
#                                  lib.size = sample(1000:10000, nrow(count.ibd), replace = T),
#                                  rel.abund  = target.rel.abund,
#                                  method.total0 = "scam")

## -----------------------------------------------------------------------------
#  simu <- Midas.sim(fitted.modified)

