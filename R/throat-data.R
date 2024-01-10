#' throat microbiome dataset
#'
#' A microbiome dataset of 60 subjects with 856 OTUs. The data were collected
#' from right and left nasopharynx and oropharynx region.
#'
#' @docType data
#'
#' @usage data(throat.otu.tab)
#'
#' @keywords datasets
#'
#' @references Charlson, E. S., Chen, J., Custers-Allen, R., Bittinger, K.,
#' Li, H., Sinha, R., Hwang, J., Bushman, F. D., & Collman, R. G. (2010).
#' Disordered microbial communities in the upper respiratory tract of cigarette smokers.
#' PloS one, 5(12), e15216. https://doi.org/10.1371/journal.pone.0015216
#'
#' @examples
#' data(throat.otu.tab)
#' \donttest{
#' MIDASim.setup(otu.tab = throat.otu.tab, mode = "nonparametric")}
"throat.otu.tab"
