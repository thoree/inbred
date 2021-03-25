#'
#' Coverage for bootstrap confidence intervals
#'
#' Currently we only consider parametric bootstrap and the kinship coefficient.
#'
#' @param ped ped object with allele frequencies.
#' @param ids Id of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param CItype Logical
#' @param conf.level Double
#' @param plot Logical
#' @param seed Integer.
#' @param verbose Logical.


#'
#' @return Returns a l
#'
#'
#' @return Returns a list with two elements, a dataframe and the result
#'   of last simulation called `lastBoot`(see documentation of `ibdBootstrap`).
#'   The columns of the data frame are
#'
#'   * `phi.hat` The estimate for each bootstrap
#'
#'   * `skew` See details
#'
#'   * `dist` See documentation of `ibdBootstrap`
#'
#'   * `lower` Lower bound for interval
#'
#'
#'   * `upper` Upper bound for interval
#'
#'   * `coverage`

#' @details See note coverage.pdf
#'
#' @export
#'
#' @examples
#' library(forrel)
#' library(ribd)
#' library(coxed) # for bca confidence intervals
#'
#' # Example Estimate coverage and more
#' ped = halfSibPed()
#' ped = setMarkers(ped, locusAttributes = NorwegianFrequencies)
#' ids = leaves(ped)
#' N = 2 # no of confidence intervals
#' B = 100 # no of  bootstraps
#' CItype = "bca"
#' conf.level = 0.8
#' plot = FALSE
#' seed = 2

#' res1 = bootPhi(ped = ped, ids = ids, N = N, B = B,
#'                CItype = CItype, conf.level = conf.level, plot = plot, seed = seed)
#'
#' boot1 = res1$lastBoot
#' phi.hat = 0.25*boot1$k1 + 0.5*boot1$k2
#' plot(density(phi.hat), main = "", xlab = "Kinship coefficient")
#' qqnorm(phi.hat)


bootPhi <- function(ped, ids = NULL, N = 2, B = 2,
                    CItype = "bca", conf.level = 0.95, plot = F, seed = NULL,
                    verbose  = F){
  if(!is.null(seed))
    set.seed(seed)
  tabParam  = tabNonparam = matrix(ncol = 6, nrow = N)
  dn = list(c(paste0("sim", 1:N)),
            c("realised", "boot", "skew", "lower", "upper", "cover"))
  dimnames(tabParam) = dimnames(tabNonparam) = dn
  boot = skew =  rep(NA, N)
  sim = profileSim(ped, N = N, ids = ids)
  kappas.obs = lapply(sim, function(x)
               as.double(ibdEstimate(x, verbose = FALSE)[1,4:6]))
  realised = unlist(lapply(kappas.obs, function(x) 0.25*x[2] + 0.5*x[3]))
  phi = kinship(ped, ids = ids)

  for (j in 1:N){
    if(verbose) cat(j, " ")
    bootPar = ibdBootstrap(ped, kappa = kappas.obs[[j]],  N = B, plot = plot)
    phis = 0.25*bootPar[,2] + 0.5*bootPar[,3]
    boot = mean(phis)
    skew = boot - realised[j]
    CI = interval(phis, method = CItype, conf.level = conf.level)
    coverage = (CI[1] <= phi) & (CI[2] >= phi)
    tabParam[j, 1:6] = c(realised[j], boot, skew, CI, coverage)

    bootNon = ibdBootstrap(sim[[j]], ids, param = "kappa", N = B,
                          method = "nonparametric", plot = plot)
    phis = 0.25*bootNon[,2] + 0.5*bootNon[,3]
    boot = mean(phis)
    skew = boot - realised[j]
    CI = interval(phis, method = CItype, conf.level = conf.level)
    coverage = (CI[1] <= phi) & (CI[2] >= phi)
    tabNonparam[j, 1:6] = c(realised[j], boot, skew, CI, coverage)
   }

  average = rbind(apply(tabParam, 2, mean),
                  apply(tabNonparam, 2, mean))
  average = data.frame(average)
  rownames(average) = c("par", "nonpar")
  list(averagedSimulations = average, pedigree.phi = phi,
       simParametric = tabParam, simNonparametric = tabNonparam)
}





