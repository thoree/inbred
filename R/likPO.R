#' Likelihood for parent offspring pair with mutation
#'
#' See Equation  and p. 25-26 and Eq. (6.14) in
#' Egeland, Kling, Mostad (2015)
#'
#' @param p Vector of allele frequencies
#' @param M Mutation matrix
#' @param AF Genotype of father
#' @param CH Genotype of child
#'
#' @details Alleles are assumed to be integers 1,2 ...
#' @examples
#'
#' ### Example p. 26 Egeland, Kling, Mostad (2015)
#' # Function:
#' a = 1; b = 2; cc = 3; d = 4 # renumbered alleles compared to p. 26
#' AF = c(a, b); CH = c(cc, d)
#' pc = 0.25; pd = 0.097826; pa = 0.25; pb = 1- pa - pd - pc
#' p = c(pa, pb, pc, pd)
#' library(pedmut) # To make mutation matrix only
#' M = mutationMatrix(model = "proportional", alleles = 1:4,
#'                            afreq = p, rate = m)
#' l1 = likPO(p, M, AF, CH)
#'
#' library(pedprobr)
#' PO = nuclearPed(father = "AF", children = "CH")
#' am = matrix(ncol = 2, c(1, 0, 3, 2, 0, 4))
#' rownames(am) = labels(PO)
#' la = list(afreq = p, alleles = 1:4, mutmod = "proportional",
#'           rate = m, name ="L1")

#' PO = setMarkers(PO,locusAttributes = la, alleleMatrix = am)
#' l2 = likelihood(PO,1)
#'
#' abs(l1-l2) < 1e-15



#' @export
likPO <- function(p,M,AF,CH){
mac <- M[AF[1],CH[1]]
mbc <- M[AF[2],CH[1]]
if(CH[1]==CH[2]){
     mad <- 0
     mbd <- 0
  } else {
     mad <- M[AF[1],CH[2]]
     mbd <- M[AF[2],CH[2]]
}
I.AF <- 1*(AF[1]==AF[2])
pAF <-I.AF*p[AF[1]]^2+2*(1-I.AF)*p[AF[1]]*p[AF[2]]
lik <- pAF*(0.5*(mac+mbc)*p[CH[2]]+0.5*(mad+mbd)*p[CH[1]])
lik
}

