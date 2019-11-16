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
#' @return A list with two elements: the likelihoods assuming independence and
#' the likelihoods for the pairs. If there is an odd number of markers,
#' the last marker is omitted for the  elements of the list.
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
#' numerator = prod(lik1[[2]])
#' Delta1 = condensedIdentity(H1, c(1, 2))
#' rho = 0.5
#' Delta2 = matrix(0, ncol = 9, nrow = 9); Delta2[9,9] = 1
#' lik2 = likPairs(a,b,cc,d, pa, pb, pc, pd, Delta = Delta1, DeltaMatrix = Delta2)
#' denominator = prod(lik2[[2]])
#' numerator/denominator
#'
#' #Alternative if pedigree is defined in `pedtools`, input can be extracted
#' ids = c(4,5)
#' nM = nMarkers(H1)
#' g = selectMarkers(H1, 1:nM)
#' odd = seq(1, nM*2, by = 2)
#' even = seq(2, nM*2,by = 2)
#' g2 = getAlleles(g)
#' a = g2[ids[1], odd]
#' b = g2[ids[1], even]
#' cc = g2[ids[2], odd]
#' d = g2[ids[2], even]
#' loci = getLocusAttributes(g)
#' pa = p[as.integer(a)] # assumes integer alleles
#' pb = p[as.integer(b)]
#' pc = p[as.integer(cc)]
#' pd = p[as.integer(d)]
#'
likPairs = function(a, b, cc, d, pa, pb, pc, pd, Delta, DeltaMatrix = NULL){
  # Checks only first argument for now,MOVE UP
  if(!is.null(DeltaMatrix) & length(a) %% 2 !=0)
    stop("Must be an even no of markers when DeltaMatrix is specified")
   nPairs = floor(length(a)/2)
   liksIndependent = likJ(a, b, cc, d, pa, pb, pc, pd, Delta)
   All = likJStates(a, b, cc, d, pa, pb, pc, pd)
   likEachPair = rep(NA, nPairs)

   if(!is.null(DeltaMatrix)){ # sequential pairs of linked markers
     for (s in 1:nPairs)
       likEachPair[s] = t(All[,2*s-1]) %*% DeltaMatrix %*% All[,2*s]
   }

   list("independentLikelihoods" = liksIndependent, "likelihoodPairs" = likEachPair)
}

