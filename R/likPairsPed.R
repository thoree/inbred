#' Likelihood for two Possibly Inbred Individuals as a function of
#' the Condensed 'Jacquard' Coefficients with `ped` input
#'
#' @param x `ped` object with markers
#' @param ids A vector oflength two identifying the individuals
#' @param Delta Double vector of length 9 summing to unity
#' @param DeltaMatrix Double 9x9 matrix of two locus identity coefficients
#' @export
#' @return The likelihood for the pairs. If there is an odd number of markers,
#' the last marker is omitted for the  elements of the list.
#' @details The implementation is based on conditioning on IBD states at two
#' linked loci. Currently, the recombination rate for each pair is assumed to be the same.
#' @seealso `ribd::condensedIdentity` and `ribd::twoLocusIdentity`.`
#' @examples
#' x = pedtools::fullSibMating(1)
#' p = c(0.4, 0.6)
#' als = 1:length(p)
#' m1 = pedtools::marker(x, alleles = als, afreq = p, "5" = 1, "6" = 1)
#' m2 = pedtools::marker(x, alleles = als, afreq = p, "5" = 1:2, "6" = 1)
#' x = pedtools::setMarkers(x, list(m1, m2))
#' rho = 0.1
#' Delta1 = ribd::condensedIdentity(x, c(5, 6))
#' Delta2 = ribd::twoLocusIdentity(x, c(5, 6), rho)
#' ids = c(5,6)
#' l1 = likPairsPed(x, ids, Delta1, Delta2)
#' l2 = pedprobr::likelihood(x, m1, m2, rho)
#' abs(l1-l2) < 1e-10

likPairsPed = function(x, ids, Delta, DeltaMatrix){
   if(!is.ped(x))
      stop("First argument must be a ped object")
   nM = nMarkers(x)
   if(nM == 0 | nM %% 2 != 0)
      stop("An even number of markers required")

   if(length(Delta) != 9  | round(sum(Delta), 6) != 1)
      stop("Second to last argument should be a vector of length 9 summing to 1")
   if(!all(dim(DeltaMatrix)==9) |  round(sum(DeltaMatrix), 6) != 1)
      stop("Last argument should a 9 by matrix with elements summing to 1")

   # Extract input for `likPairs`
   g = selectMarkers(x, 1:nM)
   odd = seq(1, nM*2, by = 2)
   even = seq(2, nM*2,by = 2)
   g2 = getAlleles(g)
   a = g2[ids[1], odd]
   b = g2[ids[1], even]
   cc = g2[ids[2], odd]
   d = g2[ids[2], even]
   loci = getLocusAttributes(g)
   pa = pb = pc = pd = rep(NA, nM)
   for (i in 1:nM){
      q = loci[[i]]$afreq
      names(q) = loci[[1]]$alleles
      pa[i] = q[a[i]]
      pb[i] = q[b[i]]
      pc[i] = q[cc[i]]
      pd[i] = q[d[i]]
   }
   likPairs(a, b, cc, d, pa, pb, pc, pd, Delta, DeltaMatrix)
}

