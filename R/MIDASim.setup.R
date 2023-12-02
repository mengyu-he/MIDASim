#' Fitting MIDAS model to microbiome data
#'
#' Midas.setup estimates parameters from a template microbiome count dataset for
#' downstream data simulation.
#'
#' @param otu.tab Numeric matrix of template microbiome count dataset. Rows are
#' samples, columns are taxa.
#' @param n.break.ties Number of replicates to break ties when ranking relative
#' abundances. Defaults to \code{100}.
#' @param mode A character indicating the modeling approach for relative abundances.
#' If \code{'parametric'}, a parametric model involving fitting a generalized
#' gamma distribution is used. If \code{'nonparametric'}, the nonparametric
#' approach involving quantile matching is applied. Note that a parametric
#' model is required if library sizes or characteristics of taxa will be modified.
#' Defaults to \code{'nonparametric'}.
#' @return Returns a list that has components:
#' \item{mat01}{Presence-absence matrix of the template data.}
#' \item{lib.size}{Observed library sizes of the template data.}
#' \item{n.taxa}{Number of taxa in the template data.}
#' \item{n.sample}{Sample size in the template data.}
#' \item{ids}{Taxa ids present in all samples in the template.}
#' \item{tetra.corr}{Estimated tetrachoric correlation of the presence-absence
#' matrix of the template.}
#' \item{corr.rel.corrected}{Estimated Pearson correlation of relative abundances,
#' transformed from Spearman's rank correlation.}
#' \item{sample.1.prop}{Proportion of non-zero cells for each subject.}
#' \item{taxa.1.prop}{Proportion of non-zeros for each taxon.}
#' \item{mean.rel.abund}{Observed mean relative abundances of each taxon.}
#' \item{rel.abund.1}{Observed non-zero relative abundances of each taxon.}
#' \item{taxa.names}{Names of taxa in the template.}
#' \item{mu.est}{Estimated location parameters for the parametric model.}
#' \item{sigma.est}{Estimated scale parameters for the parametric model.}
#' \item{Q.est}{Estimated scale parameters for the parametric model.}
#'
#' @author Mengyu He
#'
#' @examples
#'
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' # use nonparametric model
#' fitted = MIDASim.setup(otu.tab)
#'
#' # use parametric model
#' fitted = MIDASim.setup(otu.tab, mode = 'parametric')
#'
#' @export
MIDASim.setup = function(otu.tab,
                         n.break.ties = 100,
                         mode = 'nonparametric') {

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

  stopifnot("Invalid value for 'mode'. Please choose either 'parametric' or 'nonparametric'." = (mode %in% c('parametric', 'nonparametric')))

  mat01 <- ifelse(otu.tab > 0, 1, 0)
  rel.tab <- normalize_rel(otu.tab)
  obs.lib.size <- rowSums(otu.tab)
  rel.abund.1 <- apply(rel.tab, 2, function(x) x[x > 0])

  n.sample <- nrow(mat01)
  n.taxa <- ncol(mat01)
  ids <- which( colSums(mat01) == n.sample )

  sample.1.prop <- rowMeans(mat01) # proportion of 1's of samples
  taxa.1.prop <- colMeans(mat01)    # proportion of 1's of taxa
  mean.rel.abund <- colMeans(normalize_rel(otu.tab))    # mean relative abundance of taxa

  tetra.corr <- cal_tetra(mat01)

  corr.rel <- matrix(0, nrow = n.taxa, ncol = n.taxa)

  for (i in 1:n.break.ties) {
    rank.tab <- apply(rel.tab, 2,  function(x) rank(x, ties.method = "random"))
    corr.rank <- cor(rank.tab, method =  "spearman")
    corr.rel <- corr.rel + 2*sin(pi*corr.rank/6)
  }

  corr.rel <- corr.rel/n.break.ties
  corr.rel.corrected <- correct_corr(corr.rel)

  if (mode == 'parametric') {

    cv.sq = n.sample / (n.sample - 1) * (colMeans(rel.tab^2) - mean.rel.abund^2) / mean.rel.abund^2  # coefficient of variation^2
    fitted.ggamma = fit_ggamma(cv.sq, mean.rel.abund,
                               obs.lib.size, mat01)

    fitted <- append(fitted, list(mat01 = mat01,
                                  lib.size = obs.lib.size,
                                  n.taxa = n.taxa,
                                  n.sample = n.sample,
                                  ids = ids,
                                  tetra.corr = tetra.corr,
                                  sample.1.prop = sample.1.prop,
                                  taxa.1.prop = taxa.1.prop,
                                  mode = mode,
                                  corr.rel.corrected = corr.rel.corrected,
                                  cv.sq = cv.sq,
                                  mean.rel.abund = mean.rel.abund,
                                  taxa.names = colnames(otu.tab),
                                  mu.est = -fitted.ggamma$mu.est,    # mu is now (-mu)
                                  sigma.est = fitted.ggamma$sigma.est,
                                  Q.est = fitted.ggamma$Q.est
                                  ) )

  } else {

    fitted <- append(fitted, list(mat01 = mat01,
                                  lib.size = obs.lib.size,
                                  n.taxa = n.taxa,
                                  n.sample = n.sample,
                                  ids = ids,
                                  tetra.corr = tetra.corr,
                                  sample.1.prop = sample.1.prop,
                                  taxa.1.prop = taxa.1.prop,
                                  mode = mode,
                                  corr.rel.corrected = corr.rel.corrected,
                                  mean.rel.abund = mean.rel.abund,
                                  taxa.names = colnames(otu.tab),
                                  rel.abund.1 = rel.abund.1 ) )
  }

  return(fitted)

}


