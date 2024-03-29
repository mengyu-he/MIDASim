% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibd-data.R
\docType{data}
\name{count.ibd}
\alias{count.ibd}
\title{IBD microbiome dataset}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 146 rows and 614 columns.
}
\usage{
data(count.ibd)
}
\description{
A filtered microbiome dataset of patients with IBD(Inflammatory Bowel Disease)
in Human Microbiome Project 2 (HMP2).
}
\examples{
data(count.ibd)
\donttest{
MIDASim.setup(otu.tab = count.ibd, mode = "nonparametric")}
}
\references{
Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. *et al*.
Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.
*Nature* 569, 655–662 (2019). https://doi.org/10.1038/s41586-019-1237-9.
}
\keyword{datasets}
