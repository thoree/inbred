#'
#' Compare kinskip confidence intervals
#'
#' The purpose is to compare confidence intervals (CI-s) based on parametric bootstrap and
#' nonparametric bootstrap. For a simple case, iid SNPs
#' with frequency 0.5 and kappa2 = 0, the asymptotic counterparts are provided.
#'
#' @param theta Double kappa0.
#' @param x ped object with allele frequencies.
#' @param ids Id of pair.
#' @param n Integer. No of markers.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param seed Integer.


#' @return Returns a list with length 3 or 4 (if asymptotic is included).
#' The last element of the list contains averaged CI, point estimate, coverage, and dist
#' (as defined in `ibdBootstrap`),
#' the former elements contain the same information for each CI.
#'
#' @details The equations for the simple asymptotic case are documented separately.
#' CI-s for the bootstrap are calculated using the percentile method.
#' Only estimates in the asymptotic case are constrained by kappa2 = 0.
#' This makes comparison between the asymptotic and the bootstrap estimates a bit unfair.
#'
#' @export
#' @examples
#' library(forrel)
#' library(pedprobr)
#'
#' # Example 1. Assumptions for simple asymptotics met.
#' x = halfSibPed(1)
#' n = 100 # no of markers
#' N = 2 # no of simulations
#' B = 100 # no of  bootstraps
#'
#' ## Frequencies for bootstrap
#' freq = list()
#' for (i in 1:n)
#'  freq[[i]] =  list(alleles = 1:2, afreq = c(0.5, 0.5))
#' x = setMarkers(x, locusAttributes = freq)
#'
#' ## Parent offspring
#' ids = c(1, 4); theta = 0; seed = 17
#' foo1 = compareCI(theta, x, ids, n, N, B)
#'
#' ## half sibs
#' ids = c(4,5); theta = 0.5; seed = 17;
#' foo2 = compareCI(theta, x, ids, n, N, B,)
#'
#' ## Unrelated
#' ids = c(1,3); theta = 1
#' foo3 = compareCI(theta, x, ids, n, N, B)
#'
#' cbind(theta = rep(c(0,0.5,1), each = 3),
#'       rbind(foo1$average, foo2$average, foo3$average))
#'
#' # Example 2. Strong deviation from iid
#'
#' n1 = 19 # no of id markers
#' n2 = 1
#' MAF = 0.001 #minor allele frequency
#' n = n1 + n2 # no of markers
#' N = 2 # no of simulations
#' B = 100 # no of  bootstraps
#'
#' ## Frequencies for bootstrap
#' freq = list()
#' for (i in 1:n1)
#'  freq[[i]] =  list(alleles = 1:2, afreq = c(0.5, 0.5))
#' for (i in (n1+1):n)
#'  freq[[i]] =  list(alleles = 1:2, afreq = c(MAF, 1-MAF))
#' x = setMarkers(x, locusAttributes = freq)
#'
#'
#' ## Parent offspring
#' ids = c(1, 4); theta = 0; seed = 17
#' foo1 = compareCI(theta, x, ids, n, N, B, asymptotic = FALSE)
#'
#' ## half sibs
#' ids = c(4,5); theta = 0.5; seed = 17;
#' foo2 = compareCI(theta, x, ids, n, N, B, asymptotic= FALSE)
#'
#' ## Unrelated
#' ids = c(1,3); theta = 1
#' foo3 = compareCI(theta, x, ids, n, N, B, asymptotic = FALSE)
#'
#' cbind(theta = rep(c(0,0.5,1), each = 2),
#'       rbind(foo1$average, foo2$average, foo3$average))
#'
compareCI <- function(theta, x, ids, n = 2, N = 2, B = 2, seed = NULL,
                      asymptotic = TRUE){
  set.seed(seed)
  if(asymptotic){
    # Find genotype probabilities
    p = findp(x, ids)

    # Simulate from multionomial and find sufficient statistics n1, n2, n3 and n4
    tab = rmultinom(N, n, p)
    rownames(tab) = c("n11", "n12", "n13", "n21", "n22", "n23",
                      "n31", "n32", "n33" )
    n1 = tab["n11",] + tab["n33",]
    n2 = tab["n13",] + tab["n31",]
    n3 = tab["n22",]
    n4 = tab["n12",] + tab["n21",]+tab["n23",] + tab["n32",]

    ## MLE estimate and asymptotic confidence interval
    theta.MLE = pmin(c(2/(1+n1/n2)), 1)
    se = sqrt(1/(n*I1(theta.MLE)))
    CI = matrix(c("lower" = theta.MLE - 1.96*se,
                  "upper" = theta.MLE + 1.96*se), ncol = 2)
    coverage = CI[,1] <= theta & CI[,2] >= theta
    dist = sqrt(2)*abs(theta.MLE - theta)

    MLE = data.frame(lower = pmax(CI[,1],0),
                     upper = pmin(CI[,2], 1), theta.hat = theta.MLE, coverage, dist)
    MLE.average = apply(MLE, 2, mean)
  }

  CI = matrix(ncol = 2, nrow = N)
  theta.est = dist = rep(NA, N)
  # Parametric bootstrap
  for (j in 1:N){
    boot1 = ibdBootstrap(x, kappa = c(theta, 1-theta, 0),  N = B,
                         plot = F)
    CI[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.est[j] = mean(boot1$k0)
    dist[j] = mean(boot1$dist)
  }
  coverage = CI[,1] <= theta & round(CI[,2],6) >= theta
  parametric = data.frame(lower = pmax(CI[,1],0), upper = pmin(CI[,2], 1),
                          theta.hat = theta.est, coverage, dist = dist)
  parametric.average = apply(parametric,2, mean)

  # Nonparametric bootstrap
  for (j in 1:N){
    x1 = profileSim(x,  1, ids, verbose = F)[[1]]
    boot1 = ibdBootstrap(x1, ids, N = B, param = "kappa",
                        plot = F, method = "nonparametric")
    CI[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.est[j] = mean(boot1$k0)
    dist[j] = mean(boot1$dist)
  }
  coverage = CI[,1] <= theta & round(CI[,2],6) >= theta
  nonparametric = data.frame(lower = pmax(CI[,1],0), upper = pmin(CI[,2], 1),
                             theta.hat = theta.est,  coverage, dist = dist)
  nonparametric.average = apply(nonparametric,2, mean)

  # Output depends of whether asymptotic results are included
  if(asymptotic){
    average = rbind(MLE.average, parametric.average, nonparametric.average)
    row.names(average) = c("asymptotic", "parametric", "nonparametric")
    res = list(MLE = MLE,  parametric = parametric,  nonparametric = nonparametric,
               average = average)
  } else {
    average = rbind( parametric.average, nonparametric.average)
    row.names(average) = c( "parametric", "nonparametric")
    res = list(parametric = parametric,  nonparametric = nonparametric,
               average = average)
  }
  res
}


findp = function(x, ids){
  x = halfSibPed(1)
  m = marker(x, alleles = c("a", "b"), afreq = c(0.5, 0.5))
  tab0 = oneMarkerDistribution(x, ids = ids,
                               partialmarker = m, verbose = F)
  as.vector(tab0)
}

I1 = function(theta){
  p1.theta = p1(k0 = theta, k1 = 1-theta, k2 = 0)
  p2.theta = p2(k0 = theta)
  (1/256)*((1-p1.theta)/p1.theta + (1-p2.theta)/p2.theta + 2)
}

p1 = function(p = 0.5, k0 = 0.25, k1 = 0.5, k2 = 0.25)
  k0 * p^4 + k1 *p^3 + k2 * p^2

p2 = function(p = 0.5, k0 = 0.25)
  k0 * p^4

