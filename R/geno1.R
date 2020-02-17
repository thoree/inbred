#' Calculates P(G|H)
#'
#' The marker data for one individual at two markers is
#' G = (a1/a2, b1/b2). H is a n by n matrix of haplotype probabilities
#'
#' @param a1 First allele, first marker
#' @param a2 Second allele, first marker
#' @param b1 First allele, second marker
#' @param b2 Second allele, second marker

#' @export
#' @return P(G|H)
#' @examples
#' # Below I go through some cases trying to figure out
#' # what the LD probabilities should be
#'
#' p = q = c(0.1, 0.2, 0.3, 0.4)
#' # Checks when LE
#' H = p%o%q # Haplotype probs for LE
#'
#' abs(geno1(1, 1, 1, 1, H) - p[1]^2*q[1]^2) < 1e-15
#' abs(geno1(1, 2, 3, 3, H) - 2*p[1]*p[2]*q[3]^2) < 1e-15
#' abs(geno1(1, 2, 3, 4, H) - 2*p[1]*p[2]*2*q[3]*q[4]) < 1e-15
#' abs(geno1(2, 2, 4, 4, H) - p[2]^2*q[4]^2) < 1e-15
geno1 = function(a1, a2, b1, b2, H){
  I1 = a1 == a2
  I2 = b1 == b2
  H11 = H[a1, b1]
  H12 = H[a1, b2]
  H21 = H[a2, b1]
  H22 = H[a2, b2]
  pG = I1*I2*H11^2 + 2*(1-I1)*I2*H11*H21 + 2*I1*(1-I2)*H11*H12 +
      (1-I1)*(1-I2)*2*(H11*H22 + H12*H21)
  pG
}


