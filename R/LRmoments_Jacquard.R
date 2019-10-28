#' Exact expectation and variance of likelihood ratios
#'
#' Computes the expectation and variance of the likelihood ratio for two individuals,
#' assumed by HP to be related according to DeltaP and for HD to be unrelated,
#' when their true relationship is given by  DeltaT.
#' Computations apply to a single locus.
#'
#' @param DeltaP Numeric vector, Jacquard coefficients DeltaP = (Delta0,...,Delta9),
#' relationship assumed by HP.
#' @param DeltaT Numeric vector, Jacquard coefficients DeltaT = (Delta0,...,Delta9),
#' relationship assumed by HP.
#' @param p Vector of allele frequencies of a locus
#' @return Named vector with expectation and variance of LR
#' @author Hilde Kjelgaard Brustad
#' @export
#' @examples
#'
#' # Example in Section 6.4
#' DeltaP = c(0, 0, 0, 0, 0, 0, 0, 1, 0 )
#' f = 0.2
#' DeltaT = c(0, 0, f, 0, 0, 0, 0, 1-f, 0 )
#' p = c("1" = 0.2, "2" = 0.8)
#' res1 = LRmoments_Jacquard(DeltaP, DeltaT, p)
#' # Check against formulae
#' L = 2
#' mu = (L+1)/2*f + (L+3)/4*(1-f)
#' s3 = sum(1/p)
#' h883 = (3*L+s3)/4
#' h888 = (5*L+3)/8 + (s3-L)/16
#' mu2= f*h883 + (1-f)*h888
#' sigma2 = mu2 - mu^2
#' check = data.frame(ELR = c(res1[1], mu), varLR = c(res1[2], sigma2))
#' rownames(check) = c("code", "formula")
#' check

LRmoments_Jacquard <- function(DeltaP, DeltaT, p){

  # List of matrices Bi
  B = Bi_matrices(p)

  # Expected LR
  ELR = DeltaP%*%B[[9]]%*%DeltaT

  # Computes E(LR^2)
  varLR = 0
  for (i in 1:9){
    varLR =  varLR + DeltaP[i]*(DeltaP%*%B[[i]]%*%DeltaT)
  }
  # Variance of LR: varLR = E(LR^2)-E(LR)^2
  varLR = varLR-(ELR^2)

  res = c("ELR"= ELR, "VarLR" = varLR)

  return(res)
}
