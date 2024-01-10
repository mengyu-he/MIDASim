#' IBD microbiome dataset
#'
#' A filtered microbiome dataset of patients with IBD(Inflammatory Bowel Disease)
#' in Human Microbiome Project 2 (HMP2).
#'
#' @docType data
#'
#' @usage data(count.ibd)
#'
#' @keywords datasets
#'
#' @references Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. *et al*.
#' Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.
#' *Nature* 569, 655–662 (2019). https://doi.org/10.1038/s41586-019-1237-9.
#'
#' @examples
#' data(count.ibd)
#' \donttest{
#' MIDASim.setup(otu.tab = count.ibd, mode = "nonparametric")}
"count.ibd"
