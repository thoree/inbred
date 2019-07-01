#' Calculate the 9x9 matrix used to find expected value
#'
#' See inbred.pdf note
#'
#' @param p Double vector, allele frequencies
#' @export
#' @examples
#' LRH(c(0.2, 0.8))
#'
LRH <- function(p){
  L <- length(p)
  s <- sum(1/p)
  s2 <- sum(1/p^2)
  line1 <- c(s2,   s,      s,   L,      s,   L,         s,      L,   1)
  line2 <- c(s,  L^2,      L,   L,      L,   L,         L,      1,   1)
  line3 <- c(s,    L,(L+s)/2,   L,      L,   1,         L,(L+1)/2,   1)
  line4 <- c(L,    L,      L,   L,      1,   1,         1,      1,   1)
  line5 <- c(s,    L,      L,   1,(L+s)/2,   L,         L,(L+1)/2,   1)
  line6 <- c(L,    L,      1,   1,      L,   L,         1,      1,   1)
  line7 <- c(s,    L,      L,   1,      L,   1, L*(L+1)/2, (L+1)/2,  1)
  line8 <- c(L,    1,(L+1)/2,   1,(L+1)/2,   1,   (L+1)/2, (L+3)/4,  1)
  line9 <- c(1,    1,      1,   1,      1,   1,         1,       1,  1)
  rbind(line1, line2, line3, line4, line5, line6, line7, line8, line9)
}
