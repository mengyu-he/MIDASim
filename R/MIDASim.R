#' Simulating Realistic Microbiome Data using MIDAS
#'
#' Generate microbiome datasets using parameters from Midas.modify.
#'
#' @param fitted.modified Output from Midas.modify
#' @param only.rel A logical indicating whether to only simulate relative-
#' abundance data. If \code{TRUE}, then the count data will not be generated.
#' Defaults to \code{FALSE}.
#' @return Returns a list that has components:
#' \item{sim_01}{Matrix of simulated presence-absence data}
#' \item{sim_rel}{Matrix of simulated relative-abundance data}
#' \item{sim_count}{Matrix of simulated count data}
#'
#' @author Mengyu He
#'
#' @examples
#'
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' fitted = MIDASim.setup(otu.tab)
#' fitted.modified = MIDASim.modify(fitted)
#' sim = MIDASim(fitted.modified, only.rel = FALSE)
#'
#' @importFrom stats runif
#' @export
MIDASim = function(fitted.modified, only.rel = FALSE) {

  n.taxa <- fitted.modified$n.taxa
  n.sample <- fitted.modified$n.sample
  theta <- fitted.modified$theta

  ids <- union(fitted.modified$zero.id, fitted.modified$one.id)
  if (fitted.modified$mode == 'nonparametric') {
    if (length(ids) == 0) {
      mvn <- MASS::mvrnorm(n = n.sample , theta,
                           Sigma = suppressWarnings(psych::cor.smooth(fitted.modified$tetra.corr)),
                           tol = 10^-8 )
      if (n.sample == 1) mvn = matrix(mvn, nrow = 1)
      sim_01 <- ifelse( mvn >= -fitted.modified$eta, 1, 0)
    } else {
      tetra.corr <- fitted.modified$tetra.corr[-ids, -ids]
      mvn <- MASS::mvrnorm(n = n.sample , theta,
                           Sigma = suppressWarnings(psych::cor.smooth(tetra.corr)),
                           tol = 10^-8 )
      if (n.sample == 1) mvn = matrix(mvn, nrow = 1)
      sim_01 <- matrix(1, nrow = n.sample, ncol = n.taxa)
      sim_01[ ,fitted.modified$ids.left ] <- ifelse( mvn >= -fitted.modified$eta, 1, 0)
      sim_01[ ,fitted.modified$zero.id ] <- 0
    }
  } else {
    if (length(ids) == 0) {
      mvn <- MASS::mvrnorm(n = n.sample , rep(0, n.taxa),
                           Sigma = suppressWarnings(psych::cor.smooth(fitted.modified$tetra.corr)),
                           tol = 10^-8 )
      if (n.sample == 1) mvn = matrix(mvn, nrow = 1)
      mvn <- mvn - qnorm( 1 - fitted.modified$prob01.mat )
      sim_01 <- ifelse( mvn >= 0, 1, 0)
    } else {
      tetra.corr <- fitted.modified$tetra.corr[-ids, -ids]
      mvn <- MASS::mvrnorm(n = n.sample , rep(0, nrow(tetra.corr)),
                           Sigma = suppressWarnings(psych::cor.smooth(tetra.corr)),
                           tol = 10^-8 )
      if (n.sample == 1) mvn = matrix(mvn, nrow = 1)
      sim_01 <- matrix(1, nrow = n.sample, ncol = n.taxa)
      mvn <- mvn - qnorm( 1 - fitted.modified$prob01.mat[,-ids] )
      sim_01[ ,fitted.modified$ids.left ] <- ifelse( mvn >= 0, 1, 0)
      sim_01[ ,fitted.modified$zero.id ] <- 0
    }
  }

  mvn <- MASS::mvrnorm(n = n.sample, mu = rep(0, n.taxa),
                       fitted.modified$corr.rel.corrected, tol = 10^(-10))
  if (n.sample == 1) mvn = matrix(mvn, nrow = 1)

  if ( fitted.modified$mode == 'nonparametric' ) {

    rel.sim <- matrix(NA, nrow = n.sample, ncol = n.taxa)
    for (j in 1:n.taxa) {
      tmp.rank <- rep(NA, n.sample)

      ## First, rank the values generated from mvnorm for 0's
      tmp.rank[sim_01[, j]==0] <- rank(mvn[sim_01[,j]==0, j], ties.method = "random")

      ## Then, stack the ranks of values generated from mvnorm for 1's
      tmp.rank[sim_01[, j]==1] <- rank(mvn[sim_01[,j]==1, j], ties.method = "random")+sum(sim_01[, j]==0)

      n1 <- sum(sim_01[, j])        # number of 1's in simulated 0/1 from step 1

      n2 <- length(fitted.modified$rel.abund.1[[j]])            # number of 1's in the original

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

    fitted.modified$mu.est = -fitted.modified$mu.est     # mu is now (-mu)
    rel.sim = matrix(0, nrow = n.sample, ncol = n.taxa)
    for (j in 1:n.taxa) {

      if (j %in% fitted.modified$zero.id) {
        rel.sim[, j] = 0
      } else {
        params <- list(fitted.modified$mu.est[j],
                       fitted.modified$sigma.est[j],
                       fitted.modified$Q.est[j])
        p_left <- 0
        ind <- which(sim_01[,j] == 1)
        p_right <- p.gen.gamma(fitted.modified$lib.size[ind],
                               params = params)

        n1 <- sum(sim_01[, j])
        u <- runif(n1, p_left, p_right)

        rel.sim[ind, j] <- sort(q.gen.gamma( u, params = params)$pi)[rank(mvn[ind, j],
                                                                          ties.method = "random")]

      }
    }
    fitted.modified$mu.est = -fitted.modified$mu.est
  }

  sim_rel <- normalize_rel(rel.sim)

  if (only.rel) {
    colnames(sim_rel) <- colnames(sim_01) <- fitted.modified$taxa.names
    return(list(sim_01 = sim_01,
                sim_rel = sim_rel)  )
  } else {

    sim_count <- sim_rel*fitted.modified$lib.size
    sim_count <- ifelse(sim_01 == 1 & sim_count < 1, 1, round(sim_count) )

    sim_rel <- normalize_rel(sim_count)

    colnames(sim_count) <- colnames(sim_rel) <- colnames(sim_01) <- fitted.modified$taxa.names

    return(list(sim_01 = sim_01,
                sim_rel = sim_rel,
                sim_count = sim_count)  )
  }

}

