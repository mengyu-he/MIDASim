#' MOMS-PI microbiome dataset
#'
#' A filtered microbiome dataset of Multi-Omic Microbiome Study-Pregnancy
#' Initiative (MOMS-PI) in Human Microbiome Project 2 (HMP2).
#'
#' @docType data
#'
#' @usage data(count.vaginal)
#'
#' @keywords datasets
#'
#' @references Fettweis, J.M., Serrano, M.G., Brooks, J.P. et al.
#' The vaginal microbiome and preterm birth.
#' Nat Med 25, 1012–1021 (2019). https://doi.org/10.1038/s41591-019-0450-2
#'
#' @examples
#' data(count.vaginal)
#' \donttest{
#' MIDASim.setup(otu.tab = count.vaginal, mode = "nonparametric")}
"count.vaginal"
