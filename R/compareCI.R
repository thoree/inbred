#'
#' Compare bootstrap CI-s (nonparametric and parametric bootstrap, and asymptotic)
#'
#' The
#'
#' @param theta Double kappa0.
#' @param n Integer. No of markers
#' @param N Integer. No of simulations
#' @param B Integer. No of bootstraps
#' @param x ped object
#' @param ids Id of pair
#' @return Returns table.
#' @export
#' @examples
#' library(forrel)
#' library(pedprobr)
#' x = halfSibPed(1)
#' ids = c(4,5)
#' ids = c(1,3)
#' n = 100 #markers
#' N = 2 #simulations
#' B = 100 # bootstraps
#' theta = 0.5
#' theta = 1
#' seed = 17
#' foo = compareCI(theta, x, ids, n, N, B, seed)

compareCI <- function(theta, x, ids, n = 2, N = 2, B = 2, seed = NULL){
  # Find genotype probabilities
  set.seed(seed)
  p = findp(x, ids)
  # Simulate from multionomial and find sufficient statistics
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
  inside = CI[,1] <= theta & CI[,2] >= theta

  # Add point estimate and indicator for being in interval
  CI.summary.MLE = data.frame(lower = pmax(CI[,1],0),
                              upper = pmin(CI[,2], 1), theta.MLE, inside)
  coverage.MLE = sum(inside)/N

  # Frequencies for bootstrap
  freqList = list()
  for (i in 1:n)
    freqList[[i]] =  c("a" = 0.5, "b" = 0.5)

  # Parametric bootstrap
  CI.para = matrix(ncol = 2, nrow = N)
  theta.para = rep(NA, N)
  for (j in 1:N){
    boot1 = kappaBootstrap(kappa = c(theta, 1-theta, 0), N = B,
                         freqList = freqList, plot = F)
    CI.para[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.para[j] = mean(boot1[,1])
  }
  inside = CI.para[,1] <= theta & round(CI.para[,2],6) >= theta
  CI.summary.para = data.frame(lower = pmax(CI.para[,1],0),
                              upper = pmin(CI.para[,2], 1), theta.para, inside)
  coverage.para = sum(inside)/N

  # Nonparametric bootstrap
  CI.non = matrix(ncol = 2, nrow = N)
  theta.non = rep(NA, N)
  for (j in 1:N){
    boot1 = ibdBootstrap(kappa = c(theta, 1-theta, 0), N = B,
                           freqList = freqList, plot = F)
    CI.non[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.non[j] = mean(boot1[,1])
  }
  inside = CI.non[,1] <= theta & round(CI.non[,2],6) >= theta
  CI.summary.non = data.frame(lower = pmax(CI.non[,1],0),
                               upper = pmin(CI.non[,2], 1), theta.non, inside)
  coverage.non = sum(inside)/N

  list(CI.summary.MLE = CI.summary.MLE, coverage.MLE = coverage.MLE,
       CI.summary.para = CI.summary.para, coverage.para = coverage.para,
       CI.summary.non = CI.summary.non, coverage.non = coverage.non)
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

