#'
#' Find all possible genotypes for one individual
#'
#' Generates all genotypes
#'
#' @param alleles Vector of allele names
#' @return Returns matrix of all possible genotypes, one row per genotype.
#' @export

#Reworked from paramlink's allGenotypes
allGenos <- function(alleles){

  n <- length(alleles)
  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))
  v1 <- rep.int(alleles[-n],(n-1):1)
  v2 <-  NULL
  for (i in 2:n) v2 <- c(v2, alleles[i:n])
  rbind(cbind(v1=alleles, v2=alleles),cbind(v1, v2, deparse.level = 0))
}
