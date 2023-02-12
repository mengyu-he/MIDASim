#' Modifying MIDAS model
#'
#' Midas.modify modifies the fitted MIDAS.setup model according to user
#' specification that one or multiple of the following characteristics, such as
#' the library sizes, number of samples, proportions of zeros, and the taxa
#' relative abundances, can be changed. This is useful if the users wants to
#' introduce an 'effect' in simulation studies.
#'
#' @param fitted An output from Midas.setup
#' @param n.sample The number of samples to simulate. Defaults to the number of
#' samples in the template dataset.
#' @param lib.size A numeric vector of pre-specified library sizes (the length
#' should be equal to \code{n.sample} if specified).
#' @param rel.abund A numeric vector of specified mean relative abundances for
#' taxa (the length should be equal to \code{n.taxa} in \code{fitted}).
#' @param rel.abund.1 A numeric vector of specified mean relative abundances
#' among non-zero samples.
#' @param taxa.1.prop A numeric vector of specified proportion of non-zeros for
#' taxa (the length should be equal to \code{n.taxa} in \code{fitted}). Note,
#' \code{rel.abund} and \code{taxa.1.prop} can not be manipulated simultaneously.
#' @param method.total0 The method to be used in obtaining the target total
#' number of non-zeros. The default method "scam" fits a SCAM model to the
#' presence-absence data using library sizes and taxa prevalences as predictors;
#' the alternative "geometric" uses the geometric mean of the number of non-zeros
#' predicted by target library sizes and taxa prevalences.
#'
#' @return Returns a list that updates and adds the following elements to the
#' input list:
#' \item{n.sample}{The target sample size in the simulation}
#' \item{prop.1}{The updated taxa non-zero proportions}
#' \item{mu}{The mean values of the multivariate normal distribution in generating
#' presence-absence data}
#' \item{lib.size}{The target library sizes in the simulation}
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
#'                                 lib.size = sample(fitted$lib.size, nrow(otu.tab)) )
#'
#' # modify proportion of non-zeros in taxa
#' fitted.modified <- Midas.modify(fitted,
#'                                 taxa.1.prop = sample(fitted$prop.1))
#'
#' @export
Midas.modify = function(fitted,
                        lib.size = NULL,
                        rel.abund = NULL,
                        rel.abund.1 = NULL,
                        taxa.1.prop = NULL,
                        method.total0 = "scam") {

  if (!is.null(taxa.1.prop) && list(taxa.1.prop) == "same") {
    taxa.1.prop <- fitted$prop.1
  }

  if (!is.null(rel.abund) && list(rel.abund) == "same") {
    rel.abund <- fitted$obs.rel.abund
  }

  if (!is.null(rel.abund.1) && list(rel.abund.1) == "same") {
    rel.abund.1 <- unlist(lapply(fitted$non.zero.rel, mean))
  }

  arg <- c(ifelse( is.null(lib.size), 0, 1 ),
           ifelse( is.null(rel.abund), 0, 1 ),
           ifelse( is.null(rel.abund.1), 0, 1 ),
           ifelse( is.null(taxa.1.prop), 0, 1 ) )

  n.sample <- ifelse(arg[1] == 0, fitted$n.sample, length(lib.size))
  fitted$n.sample <- n.sample
  n.taxa <- fitted$n.taxa

  if ( sum(arg) == 0 ) message( "No adjustment to library size and taxa features is made." )

  stopifnot( "Only providing mean relative abundances among non-zero samples is not allowed" =  (!identical(arg[2:4], c(0, 1, 0))) )


  if ( arg[2] == 1 & arg[4] == 1 ) {
    stopifnot( "mean relative abundances must be greater or equal to proportion of non-zero cells" = (sum(rel.abund > taxa.1.prop) == 0) )
  }

  if ( arg[2] == 1 & arg[3] == 1 ) {
    stopifnot( "mean relative abundances must be smaller or equal to mean relative abundances among non-zero samples" = (sum(rel.abund > rel.abund.1) == 0) )
  }

  if ( sum(arg[2:4]) == 3 ) {
    stopifnot( "mean relative abundances must be equal to the product of mean relative abundances among non-zero samples and proportion of non-zero cells" = (identical(rel.abund, rel.abund.1 * taxa.1.prop) ))
  }

  stopifnot( "length of specified relative abundance-like quantities does not align with the original data" = all(c(length(rel.abund), length(rel.abund.1), length(taxa.1.prop)) %in% c(0, n.taxa) ) )

  Ztotal1.inv <- Ztotal2.inv <- 0
  ind <- 0

  obs.prop.1 <- fitted$prop.1
  obs.lib.size <- fitted$lib.size

  if ( sum(arg[2:4]) > 0 ) {

    if (arg[4] == 1) {

      if (arg[2] == 1) {
        rel.abund.1 <- rel.abund/fitted$prop.1
      } else if (arg[3] == 0){
        xvar <- log10(fitted$prop.1)
        scamfit.rel <- scam::scam( log(fitted$obs.rel.abund/(1-fitted$obs.rel.abund)) ~  s( xvar, bs = "mpi" ))
        fitted$rel.abund <- 1/(1+exp(-predict(scamfit.rel,
                                                  newdata = data.frame(xvar = log10(taxa.1.prop) ))))
        rel.abund.1 <- fitted$rel.abund/taxa.1.prop
        }
      fitted$prop.1 <- taxa.1.prop
    } else if (arg[3] == 1) {
      fitted$prop.1 <- rel.abund/rel.abund.1
    } else {

      if (is.null(fitted$ids)) {
        xvar <- log10(fitted$obs.rel.abund)
        yvar <- obs.prop.1
      } else {
        xvar <- log10(fitted$obs.rel.abund[-fitted$ids])
        yvar <- obs.prop.1[-fitted$ids]
      }

      scamfit.prop.1 <- scam::scam( log(yvar/(1-yvar)) ~  s( xvar, bs = "mpi" ))

      fitted$prop.1 <- 1/(1+exp(-predict(scamfit.prop.1,
                                         newdata = data.frame(xvar = log10(rel.abund) ))))

      rel.abund.1 <- rel.abund/fitted$prop.1
    }

    ind <- ind +1
    Ztotal2.inv <- 1/sum(fitted$prop.1)/n.sample

    for (j in 1:n.taxa) {
      equa <- function(x) mean(fitted$non.zero.rel[[j]]^x) - rel.abund.1[j]
      a <- uniroot(equa , c(-1000,1000))$root

      fitted$non.zero.rel[[j]] <- fitted$non.zero.rel[[j]]^a
    }
  }

  fitted$ids <- which(fitted$prop.1 == 1)
  if (length(fitted$ids) == 0) {
    fitted$mu <- qnorm(fitted$prop.1)
    ids.left <- (1:n.taxa)
  } else {
    fitted$mu <- qnorm(fitted$prop.1[-fitted$ids])
    ids.left <- (1:n.taxa)[-fitted$ids]
  }
  n.taxa2 <- n.taxa - length(fitted$ids)

  if (is.null(fitted$only.rel)) {

    if (is.null(fitted$alpha)) {
      num.1 <- obs.num.1 <- fitted$num.1

      if (arg[1] == 1) {
        xvar <- log10(obs.lib.size)
        scamfit.non0 <- scam::scam( log10(obs.num.1) ~  s( xvar, bs = "mpi" ))
        num.1 <- 10^(predict(scamfit.non0,
                             newdata = data.frame(xvar = log10(lib.size) )) )
        fitted$lib.size <- lib.size
        ind <- ind +1
        Ztotal1.inv <- 1/sum(num.1)
      }

      if (ind == 2) {
        if (method.total0 == "scam") {
          x1 <- log10(rep(obs.lib.size, n.taxa2 ))
          if (length(fitted$ids) == 0) {
            x2 <- rep(obs.prop.1 , each = n.sample )
            y <- fitted$vec01
            pred.x2 <- rep(fitted$prop.1, each = n.sample)
          } else {
            x2 <- rep(obs.prop.1[-fitted$ids] , each = n.sample )
            del <- as.vector(sapply(fitted$ids, function(x) ((x * n.sample - n.sample + 1):(x * n.sample))))
            y <- fitted$vec01[-del]
            pred.x2 <- rep(fitted$prop.1[-fitted$ids], each = n.sample)
          }

          scamfit.0 <- scam::scam(y ~ s(x1, bs = "mpi") + s(x2, bs = "mpi"), family = "binomial" )

          Ztotal <- sum(predict(scamfit.0, data.frame(x1 = log10(rep(fitted$lib.size, n.taxa2 )),
                                                      x2 = pred.x2), type = "response"))
        } else if (method.total0 == "geometric") {
          Ztotal <- 1/(Ztotal1.inv + Ztotal2.inv)
        }
      } else {
        Ztotal <- ifelse(Ztotal1.inv + Ztotal2.inv==0, sum(num.1), 1/(Ztotal1.inv + Ztotal2.inv) )
      }

      max.iter <- 100
      mu0 <- fitted$mu
      eta0 <- rep(0, n.sample)

      exeed.id <- names(which(fitted$prop.1/sum(fitted$prop.1) > n.sample/Ztotal + 10^-10))
      if (length( exeed.id) >0)
        stop("The marginal proportion of non-zeros in taxa ", paste(exeed.id, collapse = ", "), " can not be reached, please consider adjusting the quantities you specified.")

      s1 <- sum(num.1)
      s2 <- sum(fitted$prop.1)
      for (k in 1:max.iter) {
        mu <- mu0
        eta <- eta0

        eta0 <- NULL
        for (i in 1:n.sample) {
          equa <- function(x) (sum(pnorm(mu - x))+length(fitted$ids))/Ztotal - num.1[i]/s1
          eta0[i] <- pracma::fzero(fun = equa , x = c(-100,100), tol = 10^-10)$x
        }

        mu0 <- NULL
        jj <- 1
        for (j in ids.left) {
          equa <- function(x) sum(pnorm(x - eta0))/Ztotal - fitted$prop.1[j]/s2
          mu0[jj] <- pracma::fzero(fun = equa , x = c(-100,100), tol = 10^-10)$x
          jj <- jj + 1
        }

        diff <- sum(c(abs(mu0-mu), abs(eta0-eta)))
        if (diff < 10^-5 ){
          break
        }
      }
      fitted$mu <- mu0
      fitted$eta <- eta0
    } else {
      alpha <- beta <- NULL
      for (j in 1:n.taxa) {
        xbar <- unlist(lapply(fitted$non.zero.rel, mean))
        s2 <- unlist(lapply(fitted$non.zero.rel, var))

        alpha[j] <- xbar*(xbar*(1-xbar)/s2-1)
        beta[j] <- alpha[j]*(1-xbar)/xbar
      }
      fitted$alpha <- alpha
      fitted$beta <- beta
    }

  } else {
    fitted$eta <- rep(0, length( fitted$mu ) )
  }

  return(fitted)
}


