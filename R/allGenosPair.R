#'
#' Find all possible genotypes for one individual
#'
#' Generates all genotypes
#'
#' @param alleles Vector of allele names
#' @return Returns matrix of all possible genotypes, one row per genotype.
#' @export
#' @examples
#' allGenos(1:2)

#Reworked from paramlink's allGenotypes
allGenosPair <- function(alleles){
  genos1 = allGenos(alleles)
  d1 = dim(genos1)[1]
  gM = matrix(nrow = d1^2, ncol = 4)
  gM[,1:2] = matrix(rep(genos1,each = d1), ncol = 2)
  gM[,3:4] = matrix(rep(genos1,d1), ncol = 2)
  gM
}
