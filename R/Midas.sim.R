#' Simulating Realistic Microbiome Data using MIDAS
#'
#' Generate microbiome datasets using parameters from Midas.modify.
#'
#' @param fitted.modified An output from Midas.modify
#' @param only.rel A logical indicating whether to only simulate relative-
#' abundance data. If \code{TRUE}, then the count data will not be generated.
#' Defaults to \code{FALSE}.
#' @return Returns a list that has components:
#' \item{sim_01}{A matrix of simulated presence-absence data}
#' \item{sim_rel}{A matrix of simulated relative-abundance data}
#' \item{sim_count}{A matrix of simulated count data}
#'
#' @author Mengyu He
#'
#' @examples
#' library(GUniFrac)
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' fitted = Midas.setup(otu.tab)
#' fitted.modified = Midas.modify(fitted)
#' simu = Midas.sim(fitted.modified, only.rel = FALSE)
#'
#' @export
Midas.sim = function(fitted.modified, only.rel = FALSE) {

  n.taxa <- fitted.modified$n.taxa
  n.sample <- fitted.modified$n.sample
  mu <- fitted.modified$mu
  only.rel <- ifelse(is.null(fitted.modified$only.rel), only.rel, fitted.modified$only.rel)

  if (length(fitted.modified$ids) == 0) {
    mvn <- MASS::mvrnorm(n = n.sample , mu,
                         Sigma = suppressWarnings(psych::cor.smooth(fitted.modified$tetra.corr)),
                         tol = 10^-8 )
    sim_01 <- ifelse( mvn >= fitted.modified$eta, 1, 0)
  } else {
    tetra.corr <- fitted.modified$tetra.corr[-fitted.modified$ids, -fitted.modified$ids]
    mvn <- MASS::mvrnorm(n = n.sample , mu,
                         Sigma = suppressWarnings(psych::cor.smooth(tetra.corr)),
                         tol = 10^-8 )
    sim_01 <- matrix(1, nrow = n.sample, ncol = n.taxa)
    sim_01[ , -fitted.modified$ids ] <- ifelse( mvn >= fitted.modified$eta, 1, 0)
  }

  mvn <- MASS::mvrnorm(n = n.sample, mu = rep(0, n.taxa),
                       fitted.modified$corr.rel.corrected, tol = 10^(-10))

  if (is.null(fitted.modified$alpha)) {
    rel.sim <- matrix(NA, nrow = n.sample, ncol = n.taxa)
    n0 <- lengths(fitted.modified$rel.abund.1)
    for (j in 1:n.taxa) {
      tmp.rank <- rep(NA, n.sample)

      ## First, rank the values generated from mvnorm for 0's
      tmp.rank[sim_01[, j]==0] <- rank(mvn[sim_01[,j]==0, j], ties.method = "random")

      ## Then, stack the ranks of values generated from mvnorm for 1's
      tmp.rank[sim_01[, j]==1] <- rank(mvn[sim_01[,j]==1, j], ties.method = "random")+sum(sim_01[, j]==0)

      n1 <- sum(sim_01[, j])        # number of 1's in simulated 0/1 from step 1

      n2 <- n0[j]            # number of 1's in the original

      ## tmp1: non-zero relative abundances for taxon j
      tmp1 <- fitted.modified$rel.abund.1[[j]]

      ## If a taxon is assigned n1-n2 more 1's than that in the original data (n1 > n2),
      ## sample n1-n2 values from tmp1 (non-zero relative abundances) with replacement
      ## If n1 <= n2, sample n1 values from tmp1 without replacement

      if ( n1 > n2 ){
        tmp2 <- c( tmp1, sample(tmp1, n1-n2, replace = T),
                   rep(0, n.sample-n1 ))
        rel.sim[, j] <- sort(tmp2)[tmp.rank]
      } else{
        tmp2 <- c( sample(tmp1, n1, replace = FALSE), rep(0, n.sample-n1 ))
        rel.sim[, j] <- sort(tmp2)[tmp.rank]
      }
    }

  } else {
    rel.sim <- matrix(0, nrow = n.sample, ncol = n.taxa)
    for (j in 1: n.taxa) {
      tmp <- rbeta(sum(sim_01[,j]), fitted.modified$alpha[j], fitted.modified$beta[j])
      rel.sim[sim_01[,j]==1,j] <- sort(tmp)[rank(mvn[sim_01[,j]==1, j],
                                                ties.method = "random")]
    }
  }

  sim_rel <- normalize_rel(rel.sim)

  if (only.rel) {
    colnames(sim_rel) <- colnames(sim_01) <- fitted.modified$taxa.names
    return(list(sim_01 = sim_01,
                sim_rel = sim_rel)  )
  } else {

    sim_count <- sim_rel*fitted.modified$lib.size
    sim_count <- ifelse(sim_count>0 & sim_count<1, 1, round(sim_count) )

    sim_rel <- normalize_rel(sim_count)

    colnames(sim_count) <- colnames(sim_rel) <- colnames(sim_01) <- fitted.modified$taxa.names

    return(list(sim_01 = sim_01,
                sim_rel = sim_rel,
                sim_count = sim_count)  )
  }

}

