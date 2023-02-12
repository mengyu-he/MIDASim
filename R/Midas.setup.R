#' Fitting MIDAS model to microbiome data
#' Midas.setup estimates parameters from a template microbiome dataset (count data
#' or relative abundance data with zeros) for downstream data simulation.
#'
#' @param otu.tab A numeric matrix of template dataset (counts/relative
#' abundances). Rows are samples, columns are taxa.
#' @param n.break.ties The number of replicates used to break ties when ranking
#' relative abundances. Defaults to \code{100}.
#' @param extra.shrink A logical indicating whether to apply extra shrinkage
#' for the correlation matrix of presence-absence data by estimating \code{alpha}.
#' Defaults to \code{FALSE}.
#' @param fit.beta A logical indicating whether to fit a beta distribution for
#' the relative abundance of each taxon. Defaults to \code{FALSE}.
#' @return Returns a list that has components:
#' \item{vec01}{Vectorized presence-absence data}
#' \item{n.taxa}{The number of taxa in the template data.}
#' \item{ids}{Ids of taxa presenting in all samples in the template.}
#' \item{tetra.corr}{The estimated tetrachoric correlation of presence-absence
#' data of the template.}
#' \item{corr.rel.corrected}{The estimated Pearson correlation of relative
#' abundances, which is transformed from Spearman's rank correlation.}
#' \item{taxa.names}{Names of taxa in the template.}
#' \item{pbar}{The proportion of non-zeros for each subject.}
#' \item{prop.1}{The proportion of non-zeros for each taxon.}
#' \item{obs.rel.abund}{Observed mean relative abundances of each taxon.}
#' \item{non.zero.rel}{The observed non-zero relative abundances of each taxon.}
#' \item{n0}{The number of non-zeros for each taxon.}
#' \item{alpha}{If requested, the estimated \code{alpha} parameter for beta
#' distribution}
#' \item{beta}{If requested, the estimated \code{beta} parameter for beta
#' distribution}
#' \item{scamfit.lib}{The fitted Shape constrained additive model (SCAM) for
#' library size}
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
#' @export
Midas.setup = function(otu.tab,
                       n.break.ties = 100,
                       fit.beta = FALSE) {

  # otu.tab : col-taxon, row-sample
  otu.tab <- as.matrix(otu.tab)
  fitted <- list()

  if (nrow(otu.tab) > ncol(otu.tab)) {
    message("The number of samples is larger than the number of taxa, please make sure rows are samples")
  }

  if (min(rowSums(otu.tab)) == 0) {
    message("The samples with no observed taxa are removed")
    otu.tab <- otu.tab[rowSums(otu.tab)>0, ]
  }

  mat01 <- ifelse(otu.tab > 0, 1, 0)
  rel.tab <- normalize_rel(otu.tab)
  obs.lib.size <- rowSums(otu.tab)

  n.sample <- nrow(mat01)
  n.taxa <- ncol(mat01)

  if (all.equal(as.vector(obs.lib.size), rep(1, n.sample), tolerance = 10^-3 ) == TRUE) {
    message("The input matrix is a relative-abundance table.")
    only.rel <- TRUE
    fitted$only.rel <- TRUE
  }

  if ( any( colSums(mat01) == n.sample ) ) {
    ids <- which( colSums(mat01) == n.sample )
  } else ids <- NULL

  vec01 <- as.vector(mat01)
  num.1 <- rowSums(mat01) # number of 1's of samples
  prop.1 <- colMeans(mat01)    # proportion of 1's of taxa
  obs.rel.abund <- colMeans(normalize_rel(otu.tab))    # mean relative abundance of taxa

  tetra.corr <- cal_tetra(mat01)

  corr.rel <- matrix(0, nrow = n.taxa, ncol = n.taxa)

  for (i in 1:n.break.ties) {
    rank.tab <- apply(rel.tab, 2,  function(x) rank(x, ties.method = "random"))
    corr.rank <- cor(rank.tab, method =  "spearman")
    corr.rel <- corr.rel + 2*sin(pi*corr.rank/6)
  }
  corr.rel <- corr.rel/n.break.ties
  corr.rel.corrected <- correct_corr(corr.rel)

  fitted <- append(fitted, list(vec01 = vec01,
                                n.taxa = n.taxa,
                                n.sample = n.sample,
                                num.1 = num.1,
                                ids = ids,
                                tetra.corr = tetra.corr,
                                prop.1 = prop.1,
                                corr.rel.corrected = corr.rel.corrected,
                                obs.rel.abund = obs.rel.abund,
                                taxa.names = colnames(otu.tab)) )

  if (fit.beta) {
    alpha <- beta <- NULL
    for (j in 1:n.taxa) {
      tmp <- rel.tab[,j]
      tmp <- tmp[tmp>0]
      xbar <- mean(tmp)
      s2 <- var(tmp)

      alpha[j] <- xbar*(xbar*(1-xbar)/s2-1)
      beta[j] <- alpha[j]*(1-xbar)/xbar
    }

    fitted <- append(fitted, list(alpha = alpha,
                                 beta = beta) )
  }

  non.zero.rel <- list()
  n0 <- NULL
  for (j in 1:n.taxa) {
    tmp1 <- rel.tab[, j][which(rel.tab[,j]>0)]
    non.zero.rel[[j]] <- tmp1
    n0 <- c(n0, sum(mat01[,j]>0) )
  }

  fitted <- append(fitted, list(non.zero.rel = non.zero.rel,
                                n0 = n0))


  if ( is.null(fitted$only.rel) ){
    lib.size = obs.lib.size
    return( append(fitted, list(lib.size = lib.size) ))
  } else {
    return(fitted)
  }

}


