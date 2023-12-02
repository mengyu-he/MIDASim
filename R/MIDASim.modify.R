#' Modifying MIDAS model
#'
#' Midas.modify modifies the fitted MIDAS.setup model according to user
#' specification that one or multiple of the following characteristics, such as
#' the library sizes, taxa relative abundances, location parameters of the parametric
#' model can be changed. This is useful if the users wants to introduce an 'effect'
#' in simulation studies.
#'
#' @param fitted Output from Midas.setup.
#' @param lib.size Numeric vector of pre-specified library sizes (length should
#' be equal to \code{n.sample} if specified). In nonparametric mode, if \code{lib.size}
#' is specified, both \code{taxa.1.prop} and \code{sample.1.prop} should be specified.
#' @param mean.rel.abund Numeric vector of specified mean relative abundances for
#' taxa. Length should be equal to \code{n.taxa} in \code{fitted}.
#' @param gengamma.mu Numeric vector of specified location parameters for the
#' parametric model (generalized gamma model). Specify either \code{mean.rel.abund}
#' or \code{gengamma.mu}, not both. Length should be equal to \code{n.taxa} in \code{fitted}.
#' See Details. This argument is only applicable in parametric mode.
#' @param taxa.1.prop Numeric vector of specified proportion of non-zeros for
#' taxa (the length should be equal to \code{n.taxa} in \code{fitted}). This argument
#' is only applicable in nonparametric mode.
#' @param sample.1.prop Numeric vector of specified proportion of non-zeros for
#' subjects (the length should be equal to \code{n.sample} in \code{fitted}). This
#' argument is only applicable in nonparametric mode.
#'
#' @details The parametric model in MIDASim is a location-scale model, specifically, a
#' generalized gamma model for relative abundances \eqn{\pi} of a taxon. Denote \eqn{t = 1/\pi}.
#' The generalized gamma distribution for \eqn{t} is chosen so that
#'
#' \deqn{ln(t)\ =\ \mu\ +\ \sigma \cdot w}
#'
#' where \eqn{w} follows a log gamma distribution with a shape parameter \eqn{1/Q}.
#' MIDASim fits the model to the template data and estimates parameters \eqn{\mu}, \eqn{\sigma}
#' and \eqn{Q} by matching the first two moments of \eqn{\pi} and maximizing the likelihood.
#'
#' @return Returns an updated list with different elements depending on the value
#' of \code{fitted$mode}:
#' \item{n.sample}{Target sample size in the simulation.}
#' \item{lib.size}{Target library sizes in the simulation.}
#' \item{taxa.1.prop}{Updated proportions of non-zero values for each taxon.}
#' \item{sample.1.prop}{Updated proportion of non-zero cells for each subject.}
#' \item{theta}{Mean values of the multivariate normal distribution in generating
#' presence-absence data.}
#' \item{eta}{Adjustment to be applied to samples in generating presence-
#' absence data.}
#'
#' @author Mengyu He
#'
#' @examples
#'
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' fitted = MIDASim.setup(otu.tab, mode = 'parametric')
#'
#' # modify library sizes
#' fitted.modified <- MIDASim.modify(fitted,
#'                                   lib.size = sample(fitted$lib.size,
#'                                                     2*nrow(otu.tab)) )
#'
#' # modify mean relative abundances
#' fitted.modified <- MIDASim.modify(fitted,
#'                                   mean.rel.abund = fitted$mean.rel.abund * runif(fitted$n.taxa))
#' # modify location parameters of the parametric model
#' fitted.modified <- MIDASim.modify(fitted,
#'                                   gengamma.mu = fitted$mean.rel.abund * runif(fitted$n.taxa))
#'
#' @export
MIDASim.modify = function(fitted,
                          lib.size = NULL,
                          mean.rel.abund = NULL,
                          gengamma.mu = NULL,
                          sample.1.prop = NULL,
                          taxa.1.prop = NULL) {

  arg <- as.numeric( !c(is.null(lib.size), is.null(mean.rel.abund),
                        is.null(gengamma.mu), is.null(sample.1.prop),
                        is.null(taxa.1.prop)) )
  fitted$n.sample <- n.sample <- ifelse(arg[1] == 0, fitted$n.sample, length(lib.size))
  n.taxa <- fitted$n.taxa

  if ( fitted$mode == 'nonparametric' ) {
    stopifnot( "The argument gengamma.mu is not applicable in nonparametric mode " = (arg[3] == 0) )

    if (arg[1] == 1) {
      # lib.size specified
      stopifnot( "If library sizes are changed in nonparametric mode, both 'sample.1.prop' and 'taxa.1.prop' arguments must be provided. " = (sum(arg[4:5]) == 2) )
      fitted$lib.size = lib.size
    }
    if (arg[2] == 1) fitted$mean.rel.abund = mean.rel.abund
    if (arg[4] == 1) fitted$sample.1.prop = sample.1.prop
    if (arg[5] == 1) fitted$taxa.1.prop = taxa.1.prop

    stopifnot( "Ensure that 'sample.1.prop' is synchronized with 'taxa.1.prop' with respect to the total number of non-zero entries across the entire table, meaning that the product of 'sample.1.prop' and 'n.taxa' should equal the product of 'taxa.1.prop' and 'n.sample'. " = (all.equal(sum(fitted$taxa.1.prop) * fitted$n.sample, sum(fitted$sample.1.prop) * fitted$n.taxa)) )

    mean.rel.abund.1 = sapply(fitted$rel.abund.1, mean)
    if (all.equal(fitted$mean.rel.abund, mean.rel.abund.1 * fitted$taxa.1.prop) != TRUE) {

      mean.rel.abund.1 = fitted$mean.rel.abund / fitted$taxa.1.prop

      for (j in 1:n.taxa) {
        equa <- function(x) mean(fitted$rel.abund.1[[j]]^x) - mean.rel.abund.1[j]
        a <- pracma::fzero(fun = equa , x = c(-0.1, 1000))$x

        fitted$rel.abund.1[[j]] <- fitted$rel.abund.1[[j]]^a
      } # find a

    } # adjust non-zero relative abundances

    Ztotal = sum(fitted$taxa.1.prop) * n.sample

  } else {

    stopifnot( "Specify either mean.rel.abund or gengamma.mu, not both. " = (sum(arg[2:3]) < 2) )

    fitted$mu.est = -fitted$mu.est  #  mu is now (-mu)
    if (arg[2] == 1) {
      # mean rel abund specified
      fitted$mean.rel.abund = mean.rel.abund
      fitted$mu.est = mu.of.sigma.Q(sigma = fitted$sigma.est,
                                    Q = fitted$Q.est,
                                    m1 = fitted$mean.rel.abund)
    } else if (arg[3] == 1) {
      # gengamma.mu specified
      gengamma.mu = -gengamma.mu
      est = est.mu(gengamma.mu, fitted$sigma.est, fitted$Q.est)
      fitted$mu.est = gengamma.mu + log( sum(est) )
    }

    if (arg[1] == 1) {
      fitted$lib.size <- lib.size
      fitted$n.sample <- length(lib.size)
    }

    fitted$prob01.mat = get.prob01.mat(fitted$mu.est, fitted$sigma.est,
                                       fitted$Q.est, fitted$lib.size)

    sample.1 = rowSums(fitted$prob01.mat)
    taxa.1 = colSums(fitted$prob01.mat)

    fitted$sample.1.prop = sample.1/n.taxa
    fitted$taxa.1.prop = taxa.1/n.sample
    Ztotal =  sum(fitted$prob01.mat)

    fitted$mu.est = -fitted$mu.est
  }

  tmp <- check_taxa(taxa.1.prop = fitted$taxa.1.prop,
                    n.sample = n.sample,
                    Ztotal = Ztotal)
  fitted$one.id <- tmp$one.id
  fitted$zero.id <- tmp$zero.id
  Ztotal <- tmp$Ztotal

  if (length(fitted$one.id) + length(fitted$zero.id) > 0) {
    ids.left <- (1:n.taxa)[-union(fitted$one.id, fitted$zero.id)]
    n.rm <- n.taxa - length(ids.left)
    fitted$theta <- qnorm(fitted$taxa.1.prop[ids.left])
    fitted$ids.left <- ids.left
  } else {
    n.rm <- 0
    fitted$theta <- qnorm(fitted$taxa.1.prop)
    fitted$ids.left <- (1:n.taxa)
  }

  tmp = solver_mu_sigma( mu0 = fitted$theta, eta0 = rep(0, n.sample),
                         Ztotal = Ztotal,
                         sample.1.prop = fitted$sample.1.prop,
                         taxa.1.prop = fitted$taxa.1.prop,
                         ids.left = fitted$ids.left, n.sample = n.sample, n.rm = n.rm)
  fitted$theta <- tmp[["mu0"]]
  fitted$eta <- tmp[["eta0"]]

  return(fitted)
}




