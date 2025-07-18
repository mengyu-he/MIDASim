% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIDASim.setup.R
\name{MIDASim.setup}
\alias{MIDASim.setup}
\title{Fitting MIDAS model to microbiome data}
\usage{
MIDASim.setup(otu.tab, n.break.ties = 100, mode = "nonparametric")
}
\arguments{
\item{otu.tab}{Numeric matrix of template microbiome count dataset. Rows are
samples, columns are taxa.}

\item{n.break.ties}{Number of replicates to break ties when ranking relative
abundances. Defaults to \code{100}.}

\item{mode}{A character indicating the modeling approach for relative abundances.
If \code{'parametric'}, a parametric model involving fitting a generalized
gamma distribution is used. If \code{'nonparametric'}, the nonparametric
approach involving quantile matching is applied. Note that a parametric
model is required if library sizes or characteristics of taxa will be modified.
Defaults to \code{'nonparametric'}.}
}
\value{
Returns a list that has components:
\item{mat01}{Presence-absence matrix of the template data.}
\item{lib.size}{Observed library sizes of the template data.}
\item{n.taxa}{Number of taxa in the template data.}
\item{n.sample}{Sample size in the template data.}
\item{ids}{Taxa ids present in all samples in the template.}
\item{tetra.corr}{Estimated tetrachoric correlation of the presence-absence
matrix of the template.}
\item{corr.rel.corrected}{Estimated Pearson correlation of relative abundances,
transformed from Spearman's rank correlation.}
\item{sample.1.prop}{Proportion of non-zero cells for each subject.}
\item{taxa.1.prop}{Proportion of non-zeros for each taxon.}
\item{mean.rel.abund}{Observed mean relative abundances of each taxon.}
\item{rel.abund.1}{Observed non-zero relative abundances of each taxon.}
\item{taxa.names}{Names of taxa in the template.}
}
\description{
Midas.setup estimates parameters from a template microbiome count dataset for
downstream data simulation.
}
\examples{
\donttest{
data("throat.otu.tab")
  otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]

  # use nonparametric model
  fitted = MIDASim.setup(otu.tab)

  # use parametric model
  fitted = MIDASim.setup(otu.tab, mode = 'parametric')
}

}
\author{
Mengyu He
}
