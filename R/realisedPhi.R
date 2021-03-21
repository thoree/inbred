#'
#' Bootstrap intervals for realised phi
#'
#' Currently we only consider parametric bootstrap and the kinship coefficient.
#'
#' @param ped ped object with allele frequencies.
#' @param ids Id of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param seed Integer.
#' @param level Double
#' @param method Character 'parametric' or 'nonparametric'
#' @param bca Logical
#' @param plot Logical
#' @return Returns a l
#'
#' @details The equations
#'
#' @export
#' @examples
#' library(forrel)
#' library(pedprobr)
#' library(ribd)
#' library(moments)
#' library(coxed)
#'
#' # Example 1a Estimate coverage and more
#' n = 100 # no of markers
#' p = rep(0.5, n)
#' freq = list()
#' for (i in 1:n)
#'   freq[[i]] =  list(afreq = c("1" = p[i], "2" = 1- p[i]))
#' ped = quadHalfFirstCousins()
#' ids = leaves(ped)
#' ped = setMarkers(ped, locusAttributes = NorwegianFrequencies)
#' N = 1
#' B = 100
#' seed = 17
#' bca = FALSE
#' plot = T
#' date()
#' par(mfrow = c(1,2))
#' resParametric = realisedPhi(ped = ped, ids = ids, N = N, B = B,
#'                             method = "parametric", seed = seed,
#'                             level = 0.95, bca = bca, plot = plot)
#' foo = resParametric$lastBoot
#' phi.hat = 0.25*foo$k1+0.5*foo$k2
#' qqnorm(phi.hat, main = "", xlab ="")
#'
#' par(mfrow = c(1,2))
#' resNonparametric = realisedPhi(ped = ped, ids = ids, N = N, B = B,
#'                             method = "nonparametric", seed = seed, level = 0.95,
#'                             bca = bca, plot = plot)
#' foo = resNonparametric$lastBoot
#' phi.hat = 0.25*foo$k1+0.5*foo$k2
#' qqnorm(phi.hat, main = "", xlab = "")
#'
#' res = rbind(resParametric[N+1,], resNonparametric[N+1,] )
#' rownames(res) = c("parametric", "nonparametric")
#' round(res,4)
#' date()
#'
realisedPhi <- function(ped, ids = NULL, N = 2, B = 2,
                        seed = NULL, level = 0.95, method = "parametric", bca = FALSE, plot = F){
  if(!is.null(seed))
    set.seed(seed)
  n = nMarkers(ped)
  seeds = sample(1:10^9, N)
  CI = matrix(ncol = 2, nrow = N)
  phi.pedigree = kinship(ped, ids = ids)
  phi.hat = phi.realised = skew = dist = rep(NA, N)
  kappa.realised = list()
  for (j in 1:N){
      sim = markerSim(ped, N = n, verbose = F, ids = ids, seed = seeds[j])
      kappa.realised[[j]] = ibdEstimate(sim, ids = ids, param = "kappa", verbose = F)
      phi.realised[j] = 0.25*kappa.realised[[j]]$k1 + 0.5*kappa.realised[[j]]$k2
      boot1 = ibdBootstrap(sim, ids = ids, kappa = kappa.realised[[j]],
                           N = B, method = method, plot = plot)
      phis = 0.25*boot1[,2] + 0.5*boot1[,3]
      phi.hat[j] = mean(phis)
      skew[j] = skewness(phis)
      dist[j]= mean(boot1$dist)
      CI[j,] = ciRealised(phis, level = level)
  }

  coverage = (CI[,1] <= phi.pedigree) & (CI[,2] >= phi.pedigree)
  res1 = cbind(phi.realised, phi.hat, skew, dist, lower = CI[,1],
               phi.pedigree = phi.pedigree, upper = CI[,2], coverage)
  res2 = rbind(res1, apply(res1, 2, mean))
  rownames(res2) = c(paste0("sim", 1:N), "average")
  list(tab = as.data.frame(res2), lastBoot = boot1)
}

ciRealised = function(x, level = 0.95, bca = FALSE){
  if(bca)
    ci = coxed::bca(x, conf.level = level)
  else
    ci = quantile(x, c((1-level)/2,  1 - (1-level)/2))
  ci
}



