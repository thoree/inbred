#' Computes Bi, the 9 matrices in the variance formulain appendix of Brustad et al.
#'
#' @param p Vector of allele frequencies
#'
#' @return List of matrices Bi
#' @export
#' @examples
#' \dontrun{
#' Bi_matrices(c(0.2, 0.8))
#' }
Bi_matrices <- function(p){
  if(is.null(names(p)))
    names(p) = 1:length(p)
  # Generate all unordered genotypesof the pair
  genotypes = allGenosPair(names(p))
  a  = genotypes[,1]
  b  = genotypes[,2]
  cc = genotypes[,3]
  d  = genotypes[,4]
  pa = p[a]
  pb = p[b]
  pc = p[cc]
  pd = p[d]
  # Column j of below matrix to contain probabilities for genotypes given Jacquard state j
  GP = matrix(nrow = length(a), ncol = 9)
  Delta = rep(0, 9)
  for(j in 1:9){
    Delta[j] = 1
    GP[,j] = likJ(a, b, cc, d, pa, pb, pc, pd, Delta)
    Delta[j] = 0
  }
  B1 = B2 = B3 = B4 = B5 = B6 = B7 = B8 = B9 =
      matrix(nrow = 9, ncol = 9)
  for (j in 1:9){
      for (k in 1:9){
        B1[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,1])
        B2[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,2])
        B3[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,3])
        B4[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,4])
        B5[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,5])
        B6[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,6])
        B7[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,7])
        B8[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,8])
        B9[j,k] = sum(GP[,j]*GP[,k]/(GP[,9]^2)*GP[,9])
      }
    }
    res = list(B1, B2, B3, B4, B5, B6, B7, B8, B9)
  return(res)

}
