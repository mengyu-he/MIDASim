#' Modifying MIDAS model
#'
#' Midas.modify modifies the fitted MIDAS.setup model according to user
#' specification that one or multiple of the following characteristics, such as
#' the library sizes, number of samples, proportions of zeros, and the taxa
#' relative abundances, can be changed. This is useful if the users wants to
#' introduce an 'effect' in simulation studies.
#'
#' @param fitted An output from Midas.setup
#' @param lib.size A numeric vector of pre-specified library sizes (the length
#' should be equal to \code{n.sample} if specified).
#' @param mean.rel.abund A numeric vector of specified mean relative abundances
#' for taxa (the length should be equal to \code{n.taxa} in \code{fitted}).
#' @param mean.rel.abund.1 A numeric vector of specified mean relative abundances
#' among non-zero samples.
#' @param taxa.1.prop A numeric vector of specified proportion of non-zeros for
#' taxa (the length should be equal to \code{n.taxa} in \code{fitted}). If the
#' library sizes are given and \code{taxa.1.prop = 'same'}, the proportion of
#' zero cells in each taxon will remain the same as in the template data.
#'
#' @return Returns a list that updates and adds the following elements to the
#' input list:
#' \item{n.sample}{The target sample size in the simulation}
#' \item{lib.size}{The target library sizes in the simulation}
#' \item{taxa.1.prop}{The updated taxa non-zero proportions}
#' \item{sample.1.prop}{The updated proportion of non-zero cells for each subject.}
#' \item{mu}{The mean values of the multivariate normal distribution in generating
#' presence-absence data}
#' \item{eta}{The adjustment to be applied to samples in generating presence-
#' absence data}
#'
#' @author Mengyu He
#'
#' @examples
#' library(GUniFrac)
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' fitted = Midas.setup(otu.tab)
#'
#' # modify library sizes
#' fitted.modified <- Midas.modify(fitted,
#'                                 lib.size = sample(fitted$lib.size, 2*nrow(otu.tab)) )
#'
#' # modify proportion of non-zeros in taxa
#' fitted.modified <- Midas.modify(fitted,
#'                                 taxa.1.prop = sample(fitted$taxa.1.prop))
#'
#' @export
Midas.modify = function(fitted,
                        lib.size = NULL,
                        mean.rel.abund = NULL,
                        mean.rel.abund.1 = NULL,
                        taxa.1.prop = NULL) {

  arg <- c(ifelse( is.null(lib.size), 0, 1 ),
           ifelse( is.null(mean.rel.abund), 0, 1 ),
           ifelse( is.null(mean.rel.abund.1), 0, 1 ),
           ifelse( is.null(taxa.1.prop), 0, 1 ) )

  obs.n.sample <- fitted$n.sample
  fitted$n.sample <- n.sample <- ifelse(arg[1] == 0, fitted$n.sample, length(lib.size))
  n.taxa <- fitted$n.taxa
  obs.taxa.1.prop <- fitted$taxa.1.prop
  obs.lib.size <- fitted$lib.size
  obs.mean.rel.abund <- fitted$mean.rel.abund
  obs.sample.1.prop <- fitted$sample.1.prop
  obs.sample.1.ct <- obs.sample.1.prop*n.taxa

  n.same <- 0
  if (list(taxa.1.prop) == "same") {
    taxa.1.prop <- fitted$taxa.1.prop
    n.same <- n.same + 1
  }

  if (list(mean.rel.abund) == "same") {
    mean.rel.abund <- fitted$mean.rel.abund
    n.same <- n.same + 1
  }

  if (list(mean.rel.abund.1) == "same") {
    mean.rel.abund.1 <- unlist(lapply(fitted$rel.abund.1, mean))
    n.same <- n.same + 1
  }

  if (n.same >= 2 && arg[1] == 1 && obs.n.sample == n.sample) {
    # simple scaling
    fitted$lib.size <- lib.size
    arg[1] <- 0
  }

  stopifnot( "Only providing mean relative abundances among non-zero samples is not allowed" =  (!identical(arg[2:4], c(0, 1, 0))) )

  # arg[2]: mean relative abundances; arg[3]: mean relative abundances among non-zero samples; arg[4]:proportion of non-zero cells
  if ( arg[2] == 1 & arg[4] == 1 ) {
    stopifnot( "mean relative abundances must be greater or equal to proportion of non-zero cells" = (sum(mean.rel.abund > taxa.1.prop) == 0) )
  }

  if ( arg[2] == 1 & arg[3] == 1 ) {
    stopifnot( "mean relative abundances must be smaller or equal to mean relative abundances among non-zero samples" = (sum(mean.rel.abund > mean.rel.abund.1) == 0) )
  }

  if ( sum(arg[2:4]) == 3 ) {
    stopifnot( "mean relative abundances must be equal to the product of mean relative abundances among non-zero samples and proportion of non-zero cells" = (identical(mean.rel.abund, mean.rel.abund.1 * taxa.1.prop) ))
  }

  stopifnot( "length of specified relative abundance-like quantities does not align with the original data" = all(c(length(mean.rel.abund), length(mean.rel.abund.1), length(taxa.1.prop)) %in% c(0, n.taxa) ) )

  if (arg[1] == 1) {
    fitted$lib.size <- lib.size
    xvar <- log10(obs.lib.size)
    scamfit.non0 <- scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
    sample.1.ct <- 10^(predict(scamfit.non0,
                         newdata = data.frame(xvar = log10(lib.size) )) )
    sample.1.ct <- pmin(sample.1.ct, n.taxa)
    fitted$sample.1.prop <- sample.1.ct/n.taxa
  }

  if ( sum(arg[2:4]) > 0 ) {

    if (arg[4] == 1) {
      fitted$taxa.1.prop <- taxa.1.prop

      if (arg[2] == 1) {
        # full taxa information given

        fitted$mean.rel.abund <- mean.rel.abund
        mean.rel.abund.1 <- mean.rel.abund/taxa.1.prop
      } else if (arg[3] == 0){
        # only taxa.1.prop given

        xvar <- log10(obs.taxa.1.prop)
        suppressWarnings(scamfit.rel <- scam::scam( obs.mean.rel.abund ~  s( xvar, bs = "mpi" ), family = "binomial"))
        fitted$mean.rel.abund <- predict(scamfit.rel,
                                         newdata = data.frame(xvar = log10(taxa.1.prop) ),
                                         type = "response")
        mean.rel.abund.1 <- fitted$mean.rel.abund/taxa.1.prop
      }
    } else if (arg[3] == 1) {
      # full taxa information given

      fitted$taxa.1.prop <- mean.rel.abund/mean.rel.abund.1
      fitted$mean.rel.abund <- mean.rel.abund
    } else {
      # only mean.rel.abund given

      fitted$mean.rel.abund <- mean.rel.abund
      xvar <- log10(obs.mean.rel.abund)
      yvar <- obs.taxa.1.prop

      suppressWarnings(scamfit.prop.1 <- scam::scam( yvar ~  s( xvar, bs = "mpi" ), family = "binomial"))

      fitted$taxa.1.prop <- predict(scamfit.prop.1,
                                    newdata = data.frame(xvar = log10(mean.rel.abund) ),
                                    type = "response")

      mean.rel.abund.1 <- mean.rel.abund/fitted$taxa.1.prop
    }

    for (j in 1:n.taxa) {
      equa <- function(x) mean(fitted$rel.abund.1[[j]]^x) - mean.rel.abund.1[j]
      a <- uniroot(equa , c(-1000,1000))$root

      fitted$rel.abund.1[[j]] <- fitted$rel.abund.1[[j]]^a
    }
  }

  if (!is.null(fitted$alpha)) {
    tmp = lapply(fitted$rel.abund.1, alpha_beta)
    fitted <- append(fitted, list(alpha = sapply(tmp, "[", 1),
                                  beta = sapply(tmp, "[", 2)) )
    return(fitted)
  }

  fitted$ids <- which(fitted$taxa.1.prop == 1)
  n.rm <- length(fitted$ids)
  if (n.rm == 0) {
    fitted$mu <- qnorm(fitted$taxa.1.prop)
    ids.left <- (1:n.taxa)
  } else {
    fitted$mu <- qnorm(fitted$taxa.1.prop[-fitted$ids])
    ids.left <- (1:n.taxa)[-fitted$ids]
  }

  if (arg[1] == 1 && sum(arg[2:4]) == 1) {
    fitted$lib.size <- lib.size
    x1 <- log10(rep(obs.lib.size, n.taxa))
    x2 <- rep(obs.mean.rel.abund, each = obs.n.sample)
    pred.x2 <- rep(fitted$mean.rel.abund, each = n.sample)
    y <- fitted$vec01
    # joint modeling
    scamfit.0 <- scam::scam(y ~ s(x1) + s(x2),
                            family = "binomial", optimizer = "efs")

    Ztotal <- sum(predict(scamfit.0, data.frame(x1 = log10(rep(fitted$lib.size, n.taxa )),
                                                x2 = pred.x2), type = "response"))

    exeed.id <- (which(fitted$taxa.1.prop/sum(fitted$taxa.1.prop) > n.sample/Ztotal + 10^-10))
    Ztotal <- sum(sample.1.ct) - sum(fitted$taxa.1.prop[exeed.id]/sum(fitted$taxa.1.prop)*sum(sample.1.ct) - n.sample)
    if (length( exeed.id) >0) {
      fitted$ids <- union(fitted$ids, exeed.id)
      ids.left <- (1:n.taxa)[-fitted$ids]
      fitted$mu <- qnorm(fitted$taxa.1.prop[-fitted$ids])
    }
    n.rm <- n.taxa - length(ids.left)

    tmp = solver_mu_sigma( mu0 = fitted$mu, eta0 = rep(0, n.sample),
                           Ztotal = Ztotal,
                           sample.1.prop = fitted$sample.1.prop,
                           taxa.1.prop = fitted$taxa.1.prop,
                           ids.left = ids.left, n.sample = n.sample, n.rm = n.rm)
    fitted$mu <- tmp[["mu0"]]
    fitted$eta <- tmp[["eta0"]]
  } else if (arg[1] == 1 && sum(arg[2:4]) == 0) {
    # only library sizes are changed
    exeed.id <- (which(fitted$taxa.1.prop/sum(fitted$taxa.1.prop) > n.sample/sum(sample.1.ct) + 10^-10))
    Ztotal <- sum(sample.1.ct) - sum(fitted$taxa.1.prop[exeed.id]/sum(fitted$taxa.1.prop)*sum(sample.1.ct) - n.sample)
    if (length( exeed.id) >0) {
      fitted$ids <- union(fitted$ids, exeed.id)
      ids.left <- (1:n.taxa)[-fitted$ids]
      fitted$mu <- qnorm(fitted$taxa.1.prop[-fitted$ids])
    }
    n.rm <- n.taxa - length(ids.left)

    tmp = solver_mu_sigma( mu0 = fitted$mu, eta0 = rep(0, n.sample),
                           Ztotal = Ztotal,
                           sample.1.prop = fitted$sample.1.prop,
                           taxa.1.prop = obs.taxa.1.prop,
                           ids.left = ids.left, n.sample = n.sample, n.rm = n.rm)
    fitted$mu <- tmp[["mu0"]]
    fitted$eta <- tmp[["eta0"]]

  } else if (arg[1] == 1) {
    # scaling for lib.sizes, but change taxa features
    tmp = solver_mu_sigma( mu0 = fitted$mu, eta0 = rep(0, n.sample),
                           Ztotal = sum(fitted$taxa.1.prop)*n.sample,
                           sample.1.prop = fitted$sample.1.prop,
                           taxa.1.prop = fitted$taxa.1.prop,
                           ids.left = ids.left, n.sample = n.sample, n.rm = n.rm)
    fitted$mu <- tmp[["mu0"]]
    fitted$eta <- tmp[["eta0"]]
  } else {
    tmp = solver_mu_sigma( mu0 = fitted$mu, eta0 = rep(0, n.sample),
                           Ztotal = sum(fitted$taxa.1.prop)*n.sample,
                           sample.1.prop = fitted$sample.1.prop,
                           taxa.1.prop = fitted$taxa.1.prop,
                           ids.left = ids.left, n.sample = n.sample, n.rm = n.rm)
    fitted$mu <- tmp[["mu0"]]
    fitted$eta <- tmp[["eta0"]]
  }

  return(fitted)
}


