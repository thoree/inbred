#'
#' Properties of bootstrap confidence intervals for the kinship coefficient
#'
#' The purpose is to compare parametric and nonparametric confidence
#' intervals in kinship applications, currently only for the kinship coefficient.
#'
#' @param ped ped object with allele frequencies.
#' @param ids Id of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param CItype Logical. See [kinBoot::interval()]
#' @param conf.level Double
#' @param plot Logical
#' @param seed Integer.
#' @param verbose Logical.
#'
#'
#' @return A list with the following four elements:
#'
#'   * `phi`: The kinship coefficient of the pedigree.
#'
#'   * `averaged` : A data frame with two lines, one for `parametric`
#'   and one for `nonParametric`. The values are averaged over the `N` simulations.
#'   The columns `realised`, `boot`, `skew`, `lower`,  `upper`, and  `cover` are explained
#'   in `Details.`
#'
#'   * `simParametric`. The entries are as for `averaged`
#'   but for each parametric simulation.
#'
#'   * `simNonparametric` The entries are as for `averaged`
#'   but for each nonparametric simulation.

#' @details Marker data are simulated `N` times giving `N` estimates of (kappa0, kappa1, kappa2).
#' For each simulation, the `realised` phi is found. Parametric and nonparametric bootstrapping
#' is done giving `boot` (the averaged kinship coefficient from `B` bootstrap simulations)
#' and `skew` = `realised`- `boot`. There are various
#' ways to calculate the confidence interval. The default is `bca` as implemented in `coxed::bca`. The
#' variable `cover` is 1 if the confidence interval contains `phi` and 0 otherwise.
#'
#' @export
#'
#' @examples
#' # The example considers the kinship coefficients between
#' # brothers named `B1` and `B2` using the 35 markers in
#' # forrel::NorwegianFrequencies.
#'
#' library(forrel)
#' library(ribd)
#' library(coxed) # for bca confidence intervals
#'
#' ids = c("B1", "B2")
#' ped = nuclearPed(2, children = ids)
#' ped = setMarkers(ped, locusAttributes = NorwegianFrequencies)
#' N = 10 # no of confidence intervals. Increase
#' B = 100 # no of  bootstraps. Increase
#' res1 = bootPhi(ped, ids, N, B, seed = 17)
#'
#' # Basic output
#' res1[1:2]
#'
#' # Compare parametric and nonparametric estimates
#' y1 = res1$simParametric$boot
#' y2 = res1$simNonparametric$boot
#' boxplot(y1, y2, names = c("parametric", "nonparametric"),
#'         main = "Bootstrap estimates of kinship coefficient",
#'         sub = "Red stapled line: theoretical value")
#' abline(h = res1$phi, col = 'red', lty = 2)

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
  rownames(average) = c("parametric", "nonparametric")
  list(phi = phi,
       averaged = as.data.frame(average),
       simParametric = as.data.frame(tabParam),
       simNonparametric = as.data.frame(tabNonparam))
}





