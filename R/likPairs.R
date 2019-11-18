#' Likelihood for two Possibly Inbred Individuals as a function of
#' the Condensed 'Jacquard' Coefficients
#'
#' @param a Vector of positive integers (allele 1 of individual 1)
#' @param b Vector of positive integers (allele 2 of individual 1)
#' @param cc Vector of positive integers (allele 1 of individual 2)
#' @param d Vector of positive integers (allele 2 of individual 2)
#' @param pa Double vector of allele frequencies
#' @param pb Double vector of allele frequencies
#' @param pc Double vector of allele frequencies
#' @param pd Double vector of allele frequencies
#' @param Delta Double vector of length 9 summing to unity
#' @param DeltaMatrix Double 9x9 matrix of two locus identity coefficients
#' @export
#' @return The likelihoods for
#' pairs of markers. If there is an odd number of markers,
#' the last marker is omitted.
#' @details The implementation is based on conditioning on IBD states at two
#' linked loci. See `inbred::likPairsPed` for an implementation based `ped-suite` input.
#' @seealso `ribd::condensedIdentity` and `ribd::twoLocusIdentity`.`
#' @examples
#' library(pedtools)
#' library(ribd)
#' a  = c(1, 1, 1, 1)
#' b  = c(2, 2, 1, 2)
#' cc = c(1, 2, 2, 2)
#' d  = c(1, 2, 1, 1)
#' pa = p[a]
#' pb = p[b]
#' pc = p[cc]
#' pd = p[d]
#' H1 = halfSibPed(1)
#' p = c(0.4,  0.6)
#' # For plotting only
#' als = 1:length(p)
#' m = list()
#' for (i in 1:length(a))
#'    m[[i]] = marker(H1, afreq = p, alleles = als,
#'             "4" = c(a[i], b[i]), "5" = c(cc[i], d[i]) )
#' H1 = setMarkers(H1, m)
#' # plot(H1,m)
#'
#' Delta1 = condensedIdentity(H1, c(4,5))
#' rho = 0.01
#' Delta2 = twoLocusIdentity(H1, c(4,5), rho)
#' lik1 = likPairs(a,b,cc,d, pa, pb, pc, pd, Delta = Delta1, DeltaMatrix = Delta2)
#' numerator = prod(lik1)
#' Delta1 = condensedIdentity(H1, c(1, 2))
#' rho = 0.5
#' Delta2 = matrix(0, ncol = 9, nrow = 9); Delta2[9,9] = 1
#' lik2 = likPairs(a,b,cc,d, pa, pb, pc, pd, Delta = Delta1, DeltaMatrix = Delta2)
#' denominator = prod(lik2)
#' numerator/denominator

likPairs = function(a, b, cc, d, pa, pb, pc, pd, Delta, DeltaMatrix){
   nPairs = floor(length(a)/2)
   All = likJStates(a, b, cc, d, pa, pb, pc, pd)
   likEachPair = rep(NA, nPairs)
   for (s in 1:nPairs)
       likEachPair[s] = t(All[,2*s-1]) %*% DeltaMatrix %*% All[,2*s]
   likEachPair
}

