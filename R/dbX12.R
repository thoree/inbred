#' Argentinian X12 database
#'
#' The database contains frequencies for 12 X-chromosome STRs included in
#' `Investigator Argus X-12`. This is based on an Argentinian sample.
#' @name dbX12
#' @docType data
#' @keywords X12
#' @references "X-chromosome data for 12 STRs:
#' Towards an Argentinian database of forensic haplotype frequencies",
#' Garcia et al., FSI: Gen, suppl., vol. 41, pp. e8-13, 2019. DOI: https://doi.org/10.1016/j.fsigen.2019.04.005
#'
#' @examples
#' data(dbX12)
#' \dontrun{
#' library(pedtools, quietly = T)
#' x = setMarkers(singleton(1), locusAttributes = dbX12)
#' chrom(x, 1:length(db)) = "X"
#' x
#' }
NULL
