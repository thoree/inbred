#'
#' Example: bootstrap confidence intervals
#'
#' We consider parametric and nonparametric bootstrap and the kinship coefficient for some pedigrees.
#'
#' @param peds list of ped objects with allele frequencies.
#' @param idlist list of ids of of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param CItype Logical
#' @param conf.level Double
#' @param seed Integer
#' @param verbose Logical
#' @param db database
#'

#'
#' @return Returns


#' @details See note coverage.pdf
#'
#' @export
#'
#' @examples
#' peds = list(quadHalfFirstCousins(), doubleFirstCousins(), nuclearPed(2),
#'             halfSibPed(), cousinPed(1))
#' names(peds) = c("QHFC",  "DFC", "S", "H", "FC")
#' phi = unlist(lapply(peds, function(x) kinship(x, leaves(x))))
#' idlist = lapply(peds, leaves)
#' N = 1; B = 1000; seed = 1729

#' # Example 1 Many snps
#' n = 1000 # no of markers
#' p = rep(0.5, n)
#' freq = list()
#' for (i in 1:n)
#' freq[[i]] =  list(afreq = c("1" = p[i], "2" = 1- p[i]))
#' db = freq
#' res1 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db)
#'
#' # Example 2 Norwegian
#' db = NorwegianFrequencies
#' res2 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db)
#'
#' # Example Few Codis8
#' codis8 = c("CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51")
#' db = NorwegianFrequencies[codis8]
#' res3 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db)

#' res = rbind(res1, res2, res3)
#' res

examplePhi <- function(peds, idlist, N = 2, B = 10,  CItype = "bca",
                       conf.level = 0.95, seed = NULL, verbose  = TRUE, db = NULL){
  if(verbose) cat(date())
  if(!is.null(seed))
    set.seed(seed)
  n = length(peds)
  res = matrix(ncol = 6, nrow = 2*n)
  res = data.frame()
  for(i in 1:n){
    if(verbose) cat("\n", i, " ")
    peds[[i]] = setMarkers(peds[[i]], locusAttributes = db)

    res = rbind(res, bootPhi(ped = peds[[i]], ids = idlist[[i]], N = N, B = B,
                CItype = CItype, conf.level = conf.level)$average)
  }
  res = data.frame(par = rep(c("Y","N"), n), ped = rep(names(peds), each =2), res)
  rownames(res) = NULL
  if(verbose) cat(date(), "\n")
  res
}





